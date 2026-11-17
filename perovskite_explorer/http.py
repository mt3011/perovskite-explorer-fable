"""Cliente HTTP com cache em disco, repetição com espera e limite de taxa por servidor.

É a principal proteção contra interrupções: toda resposta bem-sucedida fica gravada em
SQLite assim que chega. Se a execução cair (energia, rede, Ctrl+C), ao rodar de novo
nenhuma consulta já respondida é refeita.
"""
from __future__ import annotations

import hashlib
import json
import logging
import random
import sqlite3
import time
import zlib
from pathlib import Path
from urllib.parse import urlparse

log = logging.getLogger("perovskite_explorer")

# intervalo mínimo (s) entre requisições a cada servidor
HOST_INTERVAL = {
    "pubchem.ncbi.nlm.nih.gov": 0.25,
    "www.crystallography.net": 0.5,
    "api.crossref.org": 0.25,
    "api.openalex.org": 0.2,
    "www.ebi.ac.uk": 0.25,
    "api.datacite.org": 0.25,
    "api.semanticscholar.org": 3.2,
}
_SECRET_PARAMS = {"api_key", "mailto", "email"}


class SourceError(Exception):
    """Falha de uma fonte (rede, servidor ou resposta inesperada)."""


class HttpClient:
    def __init__(self, cache_path, email: str | None = None, ttl_days: float = 30,
                 refresh: bool = False, offline: bool = False, timeout: float = 90,
                 max_tries: int = 5, session=None):
        self.cache_path = Path(cache_path)
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(str(self.cache_path), timeout=60)
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.execute("PRAGMA synchronous=NORMAL")
        self.db.execute(
            "CREATE TABLE IF NOT EXISTS cache (key TEXT PRIMARY KEY, url TEXT, status INTEGER, "
            "body BLOB, ts REAL)"
        )
        self.db.commit()
        self.email = email
        self.ttl_days = ttl_days
        self.refresh = refresh
        self.offline = offline
        self.timeout = timeout
        self.max_tries = max_tries
        self._session = session
        self._last: dict = {}
        self._fail_streak: dict = {}
        self._down: dict = {}
        self.stats = {"rede": 0, "cache": 0, "falhas": 0}

    # ------------------------------------------------------------------ cache
    @staticmethod
    def _key(method, url, params, data) -> str:
        p = sorted((k, str(v)) for k, v in (params or {}).items() if k not in _SECRET_PARAMS)
        d = sorted((k, str(v)) for k, v in data.items()) if isinstance(data, dict) else data
        raw = json.dumps([method, url, p, d], ensure_ascii=False, sort_keys=True)
        return hashlib.sha256(raw.encode()).hexdigest()

    def _cache_get(self, key, ttl_days):
        row = self.db.execute("SELECT status, body, ts FROM cache WHERE key=?", (key,)).fetchone()
        if not row:
            return None
        status, body, ts = row
        if not self.offline and ttl_days is not None and time.time() - ts > ttl_days * 86400:
            return None
        return status, zlib.decompress(body).decode("utf-8", "replace")

    def _cache_put(self, key, url, status, text):
        self.db.execute(
            "INSERT OR REPLACE INTO cache (key, url, status, body, ts) VALUES (?,?,?,?,?)",
            (key, url, status, zlib.compress(text.encode("utf-8")), time.time()),
        )
        self.db.commit()

    # ---------------------------------------------------------------- request
    @property
    def session(self):
        if self._session is None:
            import requests

            self._session = requests.Session()
            ua = "perovskite-explorer/0.1 (https://github.com/mt3011/perovskite-explorer-fable"
            ua += f"; mailto:{self.email})" if self.email else ")"
            self._session.headers["User-Agent"] = ua
        return self._session

    def request(self, url, params=None, data=None, headers=None, ttl_days=None,
                use_cache: bool = True):
        """Devolve (status, texto). Status possíveis: 200 ou 404. Demais casos -> SourceError."""
        method = "POST" if data is not None else "GET"
        key = self._key(method, url, params, data)
        ttl = self.ttl_days if ttl_days is None else ttl_days
        host = urlparse(url).netloc
        if use_cache and not self.refresh:
            hit = self._cache_get(key, ttl)
            if hit:
                self.stats["cache"] += 1
                return hit
        if self.offline:
            raise SourceError(f"modo offline e consulta ausente do cache: {url}")
        if host in self._down:
            raise SourceError(f"{host} indisponível nesta execução ({self._down[host]})")

        last_err = "?"
        for attempt in range(1, self.max_tries + 1):
            wait = HOST_INTERVAL.get(host, 0.3) - (time.monotonic() - self._last.get(host, 0))
            if wait > 0:
                time.sleep(wait)
            self._last[host] = time.monotonic()
            retry_after = None
            try:
                resp = self.session.request(method, url, params=params, data=data,
                                            headers=headers, timeout=self.timeout)
                self.stats["rede"] += 1
                status = resp.status_code
                if status == 200 or status == 404:
                    text = resp.text
                    self._fail_streak[host] = 0
                    if use_cache:
                        self._cache_put(key, url, status, text)
                    return status, text
                last_err = f"HTTP {status}"
                if status in (400, 401, 403, 405, 410, 422):
                    # erro do pedido ou de permissão: repetir não resolve
                    raise SourceError(f"{host}: {last_err} — {resp.text[:300]!r}")
                ra = resp.headers.get("Retry-After", "")
                retry_after = float(ra) if ra.replace(".", "", 1).isdigit() else None
                if status == 429 and retry_after and retry_after > 300:
                    # cota (diária) esgotada: esperar não adianta nesta execução
                    self._down[host] = "cota de uso esgotada (HTTP 429)"
                    raise SourceError(f"{host}: cota de uso esgotada (HTTP 429; libera em "
                                      f"{retry_after / 3600:.1f} h) — {resp.text[:200]!r}")
            except SourceError:
                self.stats["falhas"] += 1
                raise
            except Exception as exc:  # noqa: BLE001 - timeouts, DNS, conexão recusada, proxy...
                last_err = f"{type(exc).__name__}: {str(exc)[:200]}"
            if attempt < self.max_tries:
                delay = min(120.0, retry_after if retry_after else 2.0 ** attempt + random.random())
                log.warning("%s falhou (%s); tentativa %d/%d em %.0fs", host, last_err, attempt,
                            self.max_tries, delay)
                time.sleep(delay)

        self.stats["falhas"] += 1
        self._fail_streak[host] = self._fail_streak.get(host, 0) + 1
        if self._fail_streak[host] >= 2:
            # disjuntor: não desperdiça horas insistindo em um servidor fora do ar
            self._down[host] = last_err
            log.error("%s marcado como indisponível nesta execução (%s)", host, last_err)
        raise SourceError(f"{host}: {last_err} após {self.max_tries} tentativas")

    def get_json(self, url, params=None, data=None, headers=None, ttl_days=None):
        """JSON decodificado, ou None se o servidor respondeu 404 (não encontrado)."""
        status, text = self.request(url, params=params, data=data, headers=headers, ttl_days=ttl_days)
        if status == 404:
            return None
        text = text.strip()
        if not text:
            return None
        try:
            return json.loads(text)
        except ValueError as exc:
            # resposta corrompida não deve ficar no cache
            self.db.execute("DELETE FROM cache WHERE key=?", (self._key(
                "POST" if data is not None else "GET", url, params, data),))
            self.db.commit()
            raise SourceError(f"{urlparse(url).netloc}: resposta não é JSON ({text[:120]!r})") from exc

    def get_text(self, url, params=None, ttl_days=None):
        status, text = self.request(url, params=params, ttl_days=ttl_days)
        return None if status == 404 else text

    def close(self):
        try:
            self.db.close()
        except Exception:  # noqa: BLE001
            pass
