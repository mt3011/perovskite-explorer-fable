"""Interface de linha de comando.

  python -m perovskite_explorer buscar exemplos/organicos.txt -o resultados
  python -m perovskite_explorer relatorio -o resultados
  python -m perovskite_explorer diagnostico
  python -m perovskite_explorer analisar-cif estrutura.cif --organico "Phenazine"
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path

from . import __version__
from .chem import DEFAULT_HALOGENS, DEFAULT_METALS, EXTENDED_METALS
from .http import HttpClient, SourceError

log = logging.getLogger("perovskite_explorer")


def _setup_logging(out_dir: Path | None, verbose: bool):
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    log.addHandler(sh)
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(out_dir / "execucao.log", encoding="utf-8")
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        log.addHandler(fh)


def _config(args) -> dict:
    from .pipeline import DEFAULT_CONFIG

    cfg = dict(DEFAULT_CONFIG)
    cfg_file = getattr(args, "config", None) or ("config_local.json" if Path("config_local.json").exists() else None)
    if cfg_file:
        cfg.update(json.loads(Path(cfg_file).read_text(encoding="utf-8")))
    cfg["email"] = getattr(args, "email", None) or cfg.get("email") or os.environ.get("PEROVSKITE_EMAIL")
    cfg["openalex_api_key"] = cfg.get("openalex_api_key") or os.environ.get("OPENALEX_API_KEY")
    cfg["semanticscholar_api_key"] = cfg.get("semanticscholar_api_key") or os.environ.get("S2_API_KEY")
    if getattr(args, "metais", None):
        cfg["metais"] = [m.strip().capitalize() for m in args.metais.split(",") if m.strip()]
    if getattr(args, "metais_extras", False):
        cfg["metais"] = list(dict.fromkeys(cfg["metais"] + EXTENDED_METALS))
    if getattr(args, "halogenios", None):
        cfg["halogenios"] = [m.strip().capitalize() for m in args.halogenios.split(",") if m.strip()]
    if getattr(args, "fontes", None):
        cfg["fontes_literatura"] = [f.strip().lower() for f in args.fontes.split(",") if f.strip()]
    if getattr(args, "sem_cod", False):
        cfg["usar_cod"] = False
    if getattr(args, "sem_datacite", False):
        cfg["usar_datacite"] = False
    if getattr(args, "max_cifs", None) is not None:
        cfg["max_cifs"] = args.max_cifs
    return cfg


def _http(args, cfg) -> HttpClient:
    cache = Path(getattr(args, "cache", None) or "cache") / "http_cache.sqlite"
    return HttpClient(cache, email=cfg.get("email"), ttl_days=getattr(args, "validade_cache", 30),
                      refresh=getattr(args, "atualizar", False), offline=getattr(args, "offline", False))


def cmd_buscar(args) -> int:
    from . import pipeline

    out = Path(args.saida)
    _setup_logging(out, args.verbose)
    names = pipeline.read_list(args.lista)
    if not names:
        log.error("A lista %s está vazia.", args.lista)
        return 2
    cfg = _config(args)
    http = _http(args, cfg)
    log.info("Perovskite Explorer %s — %d orgânicos; saída em %s", __version__, len(names), out.resolve())
    if not cfg.get("email"):
        log.info("Dica: informe --email para entrar na fila 'educada' (mais rápida) do Crossref/OpenAlex.")
    try:
        results = pipeline.run(names, out, http, cfg, redo=args.refazer)
    finally:
        log.info("Requisições: %d pela rede, %d do cache, %d falhas", http.stats["rede"], http.stats["cache"],
                 http.stats["falhas"])
        http.close()
    done = len(results)
    incomplete = [r["entrada"] for r in results if not r["resumo"]["busca_completa"]]
    print(f"\nConcluído: {done}/{len(names)} orgânicos. Relatórios em: {out.resolve()}")
    print("  relatorio.html  resultados.xlsx  resumo.csv  estruturas.csv  literatura.csv  dados_completos.json")
    if done < len(names) or incomplete:
        if incomplete:
            print(f"  ATENÇÃO: busca incompleta para {len(incomplete)} orgânico(s): {', '.join(incomplete)}")
        print("  Rode o mesmo comando novamente para completar (nada do que já foi obtido é refeito).")
        return 1
    return 0


def cmd_relatorio(args) -> int:
    from . import pipeline, report

    out = Path(args.saida)
    _setup_logging(None, False)
    results = pipeline.load_results(out)
    if not results:
        log.error("Nenhum checkpoint encontrado em %s/checkpoints", out)
        return 2
    report.write_all(results, out, _config(args))
    print(f"Relatórios regenerados a partir de {len(results)} checkpoints em {out.resolve()}")
    return 0


def cmd_diagnostico(args) -> int:
    """Faz uma consulta mínima real a cada fonte e mostra o que respondeu."""
    _setup_logging(None, False)
    cfg = _config(args)
    http = _http(args, cfg)
    http.max_tries = 2
    mail = {"mailto": cfg["email"]} if cfg.get("email") else {}
    oa = dict(mail, **({"api_key": cfg["openalex_api_key"]} if cfg.get("openalex_api_key") else {}))
    checks = [
        ("PubChem", "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/pyridine/cids/JSON", None,
         lambda o: f"piridina -> CID {o['IdentifierList']['CID'][0]}"),
        ("COD", "https://www.crystallography.net/cod/result",
         {"el1": "Pb", "el2": "I", "el3": "C", "el4": "N", "strictmax": 5, "format": "json"},
         lambda o: f"{len(o)} entradas com Pb, I, C, N (até 5 elementos); campos: {', '.join(list(o[0])[:8])}…"),
        ("OpenAlex", "https://api.openalex.org/works", dict(oa, search="hybrid perovskite", **{"per-page": 1}),
         lambda o: f"{o['meta']['count']} trabalhos para 'hybrid perovskite'"),
        ("Europe PMC", "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
         {"query": '"hybrid perovskite"', "format": "json", "pageSize": 1},
         lambda o: f"{o['hitCount']} trabalhos para 'hybrid perovskite'"),
        ("Crossref", "https://api.crossref.org/works", dict(mail, **{"query.bibliographic": "hybrid perovskite", "rows": 1}),
         lambda o: f"{o['message']['total-results']} trabalhos para 'hybrid perovskite'"),
        ("DataCite/CCDC", "https://api.datacite.org/dois", {"client-id": "ccdc.csd", "page[size]": 1},
         lambda o: f"{o['meta']['total']} depósitos CCDC indexados"),
    ]
    failures = 0
    for label, url, params, show in checks:
        try:
            status, text = http.request(url, params=params, use_cache=False)
            print(f"[ OK ] {label:14s} {show(json.loads(text))}" if status == 200
                  else f"[FALHA] {label:14s} HTTP {status}")
            failures += status != 200
        except SourceError as exc:
            failures += 1
            print(f"[FALHA] {label:14s} {exc}")
        except Exception as exc:  # noqa: BLE001
            failures += 1
            print(f"[FALHA] {label:14s} resposta em formato inesperado ({type(exc).__name__}: {exc})")
    from . import csd

    print(f"[{' OK ' if csd.available() else ' -- '}] {'CSD Python API':14s} "
          + ("instalada (conector CSD será usado)" if csd.available() else "não instalada (opcional; ver guia)"))
    http.close()
    print("\nTodas as fontes responderam." if not failures else
          f"\n{failures} fonte(s) falharam. O programa funciona sem elas, mas marca a busca como incompleta.")
    return 1 if failures else 0


def cmd_analisar_cif(args) -> int:
    """Analisa um CIF local (ex.: baixado do CCDC Access Structures). Funciona offline."""
    from .cif_analysis import analyze_cif

    _setup_logging(None, False)
    query = None
    if args.organico:
        from . import pubchem

        cfg = _config(args)
        http = _http(args, cfg)
        try:
            org = pubchem.resolve(http, args.organico)
        finally:
            http.close()
        if not org.get("heavy_elements"):
            print(f"Não foi possível resolver '{args.organico}' no PubChem.")
            return 2
        query = {"heavy_elements": org["heavy_elements"], "heavy_edges": org["heavy_edges"]}
    code = 0
    for f in args.cif:
        r = analyze_cif(Path(f).read_text(encoding="utf-8", errors="replace"), query)
        print(f"\n== {f}")
        if not r.get("ok"):
            print("  não analisado:", r.get("erro"))
            code = 1
            continue
        print("  Arquitetura:", r["classificacao"])
        print(f"  Fórmula (CIF): {r['formula_cif']} | grupo espacial: {r['grupo_espacial']} | átomos na cela: {r['n_atomos_cela']}")
        print(f"  Componentes moleculares: {r['componentes']} | desordem: {'sim' if r['desordem'] else 'não'}")
        if query:
            o = r["organico"]
            print(f"  Orgânico '{args.organico}': {o['status']}" + (f" — {o.get('motivo')}" if o.get("motivo") else ""))
    return code


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="perovskite_explorer", description="Triagem de orgânicos para perovskitas híbridas "
                                "(bancos cristalográficos + literatura), sem IA.")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--cache", default="cache", help="pasta do cache HTTP (padrão: ./cache)")
        sp.add_argument("--config", help="arquivo JSON com parâmetros (padrão: config_local.json, se existir)")
        sp.add_argument("--email", help="seu e-mail, enviado ao Crossref/OpenAlex para uso da fila 'educada'")

    b = sub.add_parser("buscar", help="processa uma lista de orgânicos")
    b.add_argument("lista", help="arquivo .txt com um orgânico por linha")
    b.add_argument("-o", "--saida", default="resultados", help="pasta de saída (padrão: ./resultados)")
    b.add_argument("--metais", help=f"metais separados por vírgula (padrão: {','.join(DEFAULT_METALS)})")
    b.add_argument("--metais-extras", action="store_true", help=f"acrescenta {','.join(EXTENDED_METALS)}")
    b.add_argument("--halogenios", help=f"padrão: {','.join(DEFAULT_HALOGENS)} (acrescente F se quiser fluoretos)")
    b.add_argument("--fontes", help="fontes de literatura: openalex,europepmc,crossref,semanticscholar")
    b.add_argument("--sem-cod", action="store_true", help="não consulta o COD")
    b.add_argument("--sem-datacite", action="store_true", help="não busca depósitos CCDC via DataCite")
    b.add_argument("--max-cifs", type=int, help="máximo de CIFs verificados por orgânico (padrão: 60)")
    b.add_argument("--refazer", action="store_true", help="ignora checkpoints e reprocessa todos os orgânicos")
    b.add_argument("--atualizar", action="store_true", help="ignora o cache HTTP e consulta tudo de novo")
    b.add_argument("--offline", action="store_true", help="usa somente o cache (não acessa a rede)")
    b.add_argument("--validade-cache", type=float, default=30, help="dias de validade do cache (padrão: 30)")
    b.add_argument("-v", "--verbose", action="store_true")
    common(b)
    b.set_defaults(func=cmd_buscar)

    r = sub.add_parser("relatorio", help="regenera os relatórios a partir dos checkpoints (sem rede)")
    r.add_argument("-o", "--saida", default="resultados")
    common(r)
    r.set_defaults(func=cmd_relatorio)

    d = sub.add_parser("diagnostico", help="testa se cada fonte de dados está respondendo")
    common(d)
    d.set_defaults(func=cmd_diagnostico)

    a = sub.add_parser("analisar-cif", help="analisa arquivos CIF locais (dimensionalidade, orgânico presente?)")
    a.add_argument("cif", nargs="+")
    a.add_argument("--organico", help="nome do orgânico a procurar dentro do CIF (usa o PubChem)")
    common(a)
    a.set_defaults(func=cmd_analisar_cif)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except KeyboardInterrupt:
        print("\nInterrompido. Rode o mesmo comando para continuar de onde parou.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
