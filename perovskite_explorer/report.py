"""Relatórios: Excel (várias abas), CSVs, JSON completo e um HTML autocontido."""
from __future__ import annotations

import csv
import datetime as dt
import html
import json
from pathlib import Path

from . import __version__
from . import literature as lit
from . import structures as st


def _j(v, sep=", "):
    if v is None:
        return ""
    if isinstance(v, dict):
        return sep.join(f"{k} ({x})" for k, x in v.items())
    if isinstance(v, (list, tuple, set)):
        return sep.join(str(x) for x in v)
    return v


def _yn(v):
    return "sim" if v else "não"


def _name(r):
    return r["entrada"]


# --------------------------------------------------------------------------- #
# Tabelas (lista de cabeçalhos + linhas) — a mesma base serve a Excel, CSV e HTML
# --------------------------------------------------------------------------- #

def table_summary(results):
    head = ["Orgânico (entrada)", "Nível de exploração", "Índice de oportunidade (0–100)", "Parecer",
            "Busca completa?", "Nome no PubChem", "CID", "Fórmula", "Massa molar",
            "Estruturas híbridas confirmadas/prováveis", "…das quais confirmadas por grafo",
            "Estruturas híbridas apenas possíveis", "Compostos híbridos distintos",
            "Sais de haleto sem metal (estruturas)", "Combinações M–X reportadas", "Dimensionalidades",
            "Combinações M–X não encontradas (núcleo)", "Publicações estruturais", "Grupos independentes (aprox.)",
            "Citações somadas (estruturais)", "Primeiro ano", "Último ano", "Artigos: aditivo/passivação",
            "Artigos: ambíguos", "Artigos: menções", "Artigos não verificáveis", "Depósitos CSD ligados (DataCite)",
            "Sítios protonáveis", "Tipos de N", "Grupos ácidos", "Doadores lig. H", "Aceptores lig. H",
            "Ligações rotáveis", "XLogP", "TPSA", "Sais de haleto no PubChem", "Fator de viabilidade",
            "Termos de busca (cátion)", "Termos de busca (neutro)", "Fontes com falha", "Avisos"]
    rows = []
    for r in results:
        o, s = r.get("organico") or {}, r["resumo"]
        salts = (o.get("formas_relacionadas") or {}).get("sais_de_haleto") or []
        rows.append([
            _name(r), s["nivel_exploracao"], s["indice_oportunidade"], s["parecer"], _yn(s["busca_completa"]),
            o.get("titulo", ""), o.get("cid", ""), o.get("formula_hill", ""), o.get("massa_molar", ""),
            s["n_estruturas_hibridas_fortes"], s["n_estruturas_hibridas_confirmadas"],
            s["n_estruturas_hibridas_possiveis"], s["n_compostos_hibridos"], s["n_sais_haleto_sem_metal"],
            _j(s["combinacoes_MX"]), _j(s["dimensionalidades"]), _j(s["lacunas_MX"]),
            s["n_publicacoes_estruturais"], s["n_grupos_independentes"], s["citacoes_estruturais"],
            s["primeiro_ano"] or "", s["ultimo_ano"] or "", s["n_artigos_aditivo"], s["n_artigos_ambiguos"],
            s["n_mencoes"], s["n_nao_verificaveis"], s["n_depositos_csd_ligados"],
            o.get("n_sitios_protonaveis", ""),
            _j([f"{x['tipo']} [{x['basicidade']}]" for x in o.get("sitios_N", [])], "; "),
            _j(sorted(set(o.get("grupos_acidos", [])))), o.get("doadores_lig_H", ""), o.get("aceptores_lig_H", ""),
            o.get("ligacoes_rotaveis", ""), o.get("xlogp", ""), o.get("tpsa", ""),
            _j([f"{x['nome']} (CID {x['cid']})" for x in salts[:6]], "; "), s["fator_viabilidade"],
            _j((r.get("termos") or {}).get("cations", []), "; "), _j((r.get("termos") or {}).get("neutros", []), "; "),
            _j(s["fontes_com_falha"]), _j(o.get("avisos", []), " | "),
        ])
    rows.sort(key=lambda x: (-x[2], x[0].lower()))
    return head, rows


def table_compounds(results):
    head = ["Orgânico", "Composto (fórmula reduzida, sem H)", "Categoria", "Melhor confiança", "Metais", "Haletos",
            "Nº de depósitos", "Nº de publicações (DOIs)", "Primeiro ano", "Último ano", "Grupos espaciais",
            "Arquitetura estrutural", "Entradas", "DOIs"]
    rows = []
    for r in results:
        for g in r.get("compostos", []):
            rows.append([_name(r), g["composto"], g["categoria"], g["melhor_confianca"], _j(g["metais"]),
                         _j(g["haletos"]), g["n_depositos"], g["n_publicacoes"], g["primeiro_ano"] or "",
                         g["ultimo_ano"] or "", _j(g["grupos_espaciais"]), _j(g["classificacoes"], " | "),
                         _j(g["entradas"]), _j(g["dois"])])
    return head, rows


def table_structures(results):
    head = ["Orgânico", "Fonte", "ID", "Link", "Confiança", "Motivo da confiança", "Categoria", "Fórmula",
            "Nome químico no banco", "Arquitetura estrutural", "Dimensionalidade M–X", "Compartilhamento",
            "Nº coordenação M–X", "d(M–X) média (Å)", "Orgânico coordenado ao metal?", "Desordem no CIF?",
            "Orgânicos por fórmula (k)", "Estado do orgânico (balanço de H)", "Espécies extras", "Grupo espacial",
            "a", "b", "c", "α", "β", "γ", "Volume", "Z", "T (K)", "R", "Ano", "Revista", "DOI", "Título do artigo",
            "Encontrada por"]
    rows = []
    for r in results:
        for s in r.get("estruturas", []):
            a = s.get("analise") or {}
            m = s.get("casamento_formula") or {}
            ok = a.get("ok")
            dim = a.get("dimensionalidade") if ok else None
            cell = (s.get("cela") or [None] * 6)
            if not any(cell) and ok:
                cell = a.get("cela") or cell
            rows.append([
                _name(r), s["fonte"], s["id"], s.get("url", ""), s["confianca"], s.get("motivo_confianca", ""),
                s["categoria"], s.get("formula_txt", ""), s.get("nome", ""), s.get("classificacao", ""),
                f"{dim}D" if dim is not None else "", _j(a.get("compartilhamento") if ok else ""),
                _j(a.get("numeros_de_coordenacao_MX") if ok else ""), a.get("distancia_MX_media", "") if ok else "",
                _yn(a.get("organico_coordenado_ao_metal")) if ok else "", _yn(a.get("desordem")) if ok else "",
                m.get("k", ""), m.get("papel", ""), _j([f"{x['n']}×{x['especie']}" for x in m.get("extras", [])]),
                s.get("grupo_espacial", ""), *[c if c is not None else "" for c in cell[:6]],
                s.get("volume") or "", s.get("Z") or "", s.get("temperatura_K") or "", s.get("R") or "",
                s.get("ano") or "", s.get("revista", ""), s.get("doi") or "", s.get("titulo", ""),
                _j(s.get("origem", [])),
            ])
    return head, rows


def table_literature(results):
    head = ["Orgânico", "Uso do orgânico", "Relevante?", "Ano", "Título", "Revista", "Autores", "Citações", "DOI",
            "Link", "Fontes", "Termos no título", "Termos no resumo", "Pistas estruturais", "Pistas de aditivo",
            "Pontos estrutural", "Pontos aditivo", "Trecho que cita o orgânico", "Estruturas ligadas",
            "Depósitos CCDC do artigo", "Tem resumo?", "Tipo"]
    rows = []
    for r in results:
        for p in r.get("literatura", []):
            authors = p.get("autores") or []
            rows.append([
                _name(r), p.get("uso", ""), _yn(p.get("relevante")), p.get("ano") or "", p.get("titulo", ""),
                p.get("revista", ""), _j(authors[:6]) + (" et al." if len(authors) > 6 else ""),
                p.get("citacoes") if p.get("citacoes") is not None else "", p.get("doi") or "", p.get("url", ""),
                _j(p.get("fontes", [])), _j(p.get("termos_no_titulo", [])), _j(p.get("termos_no_resumo", [])),
                _j(p.get("pistas_estruturais", []), "; "), _j(p.get("pistas_aditivo", []), "; "),
                p.get("pontos_estrutural", ""), p.get("pontos_aditivo", ""), p.get("trecho", ""),
                _j(p.get("estruturas_ligadas", [])), _j(p.get("depositos_csd", [])), _yn(p.get("tem_resumo")),
                p.get("tipo", ""),
            ])
    return head, rows


def table_csd(results):
    head = ["Orgânico", "Nº CCDC", "DOI do depósito", "DOI do artigo", "Título do depósito", "Ano", "Detalhes"]
    rows = [[_name(r), d.get("ccdc") or "", d.get("doi_deposito") or "", d.get("doi_artigo") or "",
             d.get("titulo", ""), d.get("ano") or "", d.get("detalhes", "")]
            for r in results for d in r.get("depositos_csd", [])]
    return head, rows


def table_discarded(results):
    head = ["Orgânico", "Fonte", "ID", "Link", "Fórmula", "Nome químico", "Motivo do descarte", "DOI"]
    rows = [[_name(r), s["fonte"], s["id"], s.get("url", ""), s.get("formula_txt", ""), s.get("nome", ""),
             s.get("motivo_confianca", ""), s.get("doi") or ""]
            for r in results for s in r.get("descartadas", [])]
    return head, rows


def table_sources(results):
    head = ["Orgânico", "Fonte", "Funcionou?", "Detalhes"]
    rows = []
    for r in results:
        for k, v in (r.get("fontes") or {}).items():
            det = "; ".join(f"{a}={_j(b) if not isinstance(b, (int, float, str)) else b}" for a, b in v.items() if a != "ok")
            rows.append([_name(r), k, _yn(v.get("ok")), det])
    return head, rows


def table_params(cfg, results, total):
    head = ["Parâmetro", "Valor"]
    hidden = {"openalex_api_key", "semanticscholar_api_key", "email"}
    rows = [["Versão do programa", __version__],
            ["Gerado em", dt.datetime.now().strftime("%Y-%m-%d %H:%M")],
            ["Orgânicos na lista", total if total is not None else len(results)],
            ["Orgânicos processados", len(results)]]
    rows += [[k, _j(v)] for k, v in cfg.items() if k not in hidden]
    return head, rows


SHEETS = [
    ("Resumo", table_summary), ("Compostos", table_compounds), ("Estruturas", table_structures),
    ("Literatura", table_literature), ("Depósitos CSD", table_csd), ("Descartadas", table_discarded),
    ("Fontes", table_sources),
]


# --------------------------------------------------------------------------- #
# Escrita
# --------------------------------------------------------------------------- #

def write_all(results, out_dir, cfg, total=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tables = {name: fn(results) for name, fn in SHEETS}
    tables["Parâmetros"] = table_params(cfg, results, total)
    _write_json(results, out_dir / "dados_completos.json")
    for name, fname in (("Resumo", "resumo.csv"), ("Compostos", "compostos.csv"),
                        ("Estruturas", "estruturas.csv"), ("Literatura", "literatura.csv")):
        _write_csv(tables[name], out_dir / fname)
    _write_xlsx(tables, out_dir / "resultados.xlsx")
    _replace(out_dir / "relatorio.html", render_html(results, cfg, total))


def _replace(path: Path, text: str):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def _write_json(results, path: Path):
    from .pipeline import _default

    _replace(path, json.dumps(results, ensure_ascii=False, indent=1, default=_default))


def _write_csv(table, path: Path):
    head, rows = table
    tmp = path.with_suffix(".csv.tmp")
    with open(tmp, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh, delimiter=";")
        w.writerow(head)
        w.writerows(rows)
    tmp.replace(path)


def _write_xlsx(tables, path: Path):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    wb.remove(wb.active)
    fills = {"INEXPLORADO": "C6EFCE", "SEM USO": "DDEBF7", "INDÍCIOS": "FFF2CC", "POUCO": "FFE699",
             "MODERADAMENTE": "F8CBAD", "BEM": "D9D9D9"}
    for name, (head, rows) in tables.items():
        ws = wb.create_sheet(name[:31])
        ws.append(head)
        for row in rows:
            ws.append([_cell(v) for v in row])
        for c in ws[1]:
            c.font = Font(bold=True, color="FFFFFF")
            c.fill = PatternFill("solid", fgColor="1F3864")
            c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.freeze_panes = "B2"
        if rows:
            ws.auto_filter.ref = ws.dimensions
        for i, h in enumerate(head, 1):
            longest = max([len(str(h))] + [len(str(r[i - 1])) for r in rows[:200]])
            ws.column_dimensions[get_column_letter(i)].width = max(10, min(60, longest + 2))
        if name == "Resumo":
            for row in ws.iter_rows(min_row=2):
                level = str(row[1].value or "")
                for key, color in fills.items():
                    if level.startswith(key):
                        row[1].fill = PatternFill("solid", fgColor=color)
                row[3].alignment = Alignment(wrap_text=True, vertical="top")
            ws.column_dimensions["D"].width = 90
        link_cols = [i for i, h in enumerate(head) if h == "Link"]
        for row in ws.iter_rows(min_row=2):
            for i in link_cols:
                if str(row[i].value or "").startswith("http"):
                    row[i].hyperlink = row[i].value
                    row[i].font = Font(color="0563C1", underline="single")
    tmp = path.with_suffix(".xlsx.tmp")
    wb.save(tmp)
    tmp.replace(path)


def _cell(v):
    if isinstance(v, str):
        # caracteres de controle quebram o XML do Excel
        v = "".join(ch for ch in v if ch >= " " or ch in "\t\n")
        return v[:32000]
    return v


# --------------------------------------------------------------------------- #
# HTML
# --------------------------------------------------------------------------- #

_CSS = """
:root{--bg:#fff;--fg:#1c2330;--muted:#5b6575;--line:#d9dee7;--card:#f6f8fb;--accent:#1f4e9c;
--g:#1e7b3a;--b:#1f4e9c;--y:#946200;--o:#b3541e;--gr:#5b6575}
@media(prefers-color-scheme:dark){:root{--bg:#12161d;--fg:#e6e9ef;--muted:#9aa4b4;--line:#2b3340;--card:#1a2029;
--accent:#8ab4f8;--g:#6fcf8a;--b:#8ab4f8;--y:#e0b354;--o:#f0986a;--gr:#9aa4b4}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,Segoe UI,Roboto,sans-serif}
main{max-width:1180px;margin:0 auto;padding:24px 16px 64px}h1{font-size:1.6rem;margin:.2em 0}h2{font-size:1.2rem;margin:1.6em 0 .5em}
h3{font-size:1rem;margin:1.2em 0 .4em}.muted{color:var(--muted)}a{color:var(--accent)}
table{border-collapse:collapse;width:100%;font-size:.88rem}th,td{border-bottom:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}
th{font-weight:600;white-space:nowrap}.scroll{overflow-x:auto}.num{text-align:right;font-variant-numeric:tabular-nums}
.tag{display:inline-block;padding:1px 8px;border-radius:10px;border:1px solid currentColor;font-size:.78rem;font-weight:600;white-space:nowrap}
.INEXPLORADO{color:var(--g)}.SEM{color:var(--b)}.INDICIOS{color:var(--y)}.POUCO{color:var(--y)}.MODERADAMENTE{color:var(--o)}.BEM{color:var(--gr)}
details{background:var(--card);border:1px solid var(--line);border-radius:8px;margin:10px 0;padding:0 14px}
details[open]{padding-bottom:12px}summary{cursor:pointer;padding:10px 0;font-weight:600}
.warn{color:var(--o);font-weight:600}.kv{display:grid;grid-template-columns:max-content 1fr;gap:2px 14px;font-size:.9rem}
.kv dt{color:var(--muted)}.kv dd{margin:0}.small{font-size:.82rem}
"""


def _e(v):
    return html.escape(str(v if v is not None else ""))


def _tag(level):
    cls = level.split()[0].replace("Í", "I")
    return f'<span class="tag {cls}">{_e(level)}</span>'


def _doi_link(doi):
    return f'<a href="https://doi.org/{_e(doi)}">{_e(doi)}</a>' if doi else ""


def render_html(results, cfg, total=None) -> str:
    total = total if total is not None else len(results)
    ordered = sorted(results, key=lambda r: (-r["resumo"]["indice_oportunidade"], r["entrada"].lower()))
    incomplete = [r["entrada"] for r in results if not r["resumo"]["busca_completa"]]
    parts = ["<!doctype html><html lang='pt-BR'><head><meta charset='utf-8'>",
             "<meta name='viewport' content='width=device-width,initial-scale=1'>",
             "<title>Perovskite Explorer — relatório</title>", f"<style>{_CSS}</style></head><body><main>",
             "<h1>Perovskite Explorer — triagem de orgânicos</h1>",
             f"<p class='muted'>Gerado em {dt.datetime.now():%d/%m/%Y %H:%M} · versão {__version__} · "
             f"{len(results)} de {total} orgânicos processados · metais: {_e(_j(cfg.get('metais', [])))} · "
             f"haletos: {_e(_j(cfg.get('halogenios', [])))}</p>"]
    if len(results) < total:
        parts.append(f"<p class='warn'>Execução parcial: faltam {total - len(results)} orgânicos. "
                     "Rode o mesmo comando para continuar.</p>")
    if incomplete:
        parts.append(f"<p class='warn'>Busca incompleta (alguma fonte falhou) para: {_e(', '.join(incomplete))}. "
                     "Não trate esses como inexplorados antes de rodar novamente.</p>")
    parts.append("<h2>Resumo (ordenado pelo índice de oportunidade)</h2><div class='scroll'><table><tr>"
                 "<th>Orgânico</th><th>Fórmula</th><th>Nível de exploração</th><th class='num'>Índice</th>"
                 "<th class='num'>Estruturas</th><th class='num'>Publ. estruturais</th><th class='num'>Grupos</th>"
                 "<th class='num'>Aditivo</th><th>M–X reportados</th><th>Dimens.</th><th class='num'>Sítios básicos</th></tr>")
    for n, r in enumerate(ordered):
        o, s = r.get("organico") or {}, r["resumo"]
        parts.append(
            f"<tr><td><a href='#o{n}'>{_e(r['entrada'])}</a>{'' if s['busca_completa'] else ' <span class=warn>⚠</span>'}</td>"
            f"<td>{_e(o.get('formula_hill', ''))}</td><td>{_tag(s['nivel_exploracao'])}</td>"
            f"<td class='num'>{s['indice_oportunidade']}</td><td class='num'>{s['n_estruturas_hibridas_fortes']}"
            f"{'+' + str(s['n_estruturas_hibridas_possiveis']) + '?' if s['n_estruturas_hibridas_possiveis'] else ''}</td>"
            f"<td class='num'>{s['n_publicacoes_estruturais']}</td><td class='num'>{s['n_grupos_independentes']}</td>"
            f"<td class='num'>{s['n_artigos_aditivo']}</td><td>{_e(_j(s['combinacoes_MX']))}</td>"
            f"<td>{_e(_j(s['dimensionalidades']))}</td><td class='num'>{_e(o.get('n_sitios_protonaveis', ''))}</td></tr>")
    parts.append("</table></div><p class='small muted'>Estruturas: confirmadas/prováveis (+N? = apenas possíveis). "
                 "Todos os detalhes e colunas adicionais estão em <b>resultados.xlsx</b>.</p><h2>Detalhes por orgânico</h2>")
    for n, r in enumerate(ordered):
        parts.append(_organic_block(n, r))
    parts.append("</main></body></html>")
    return "".join(parts)


def _organic_block(n, r) -> str:
    o, s = r.get("organico") or {}, r["resumo"]
    p = [f"<details id='o{n}'><summary>{_e(r['entrada'])} — {_tag(s['nivel_exploracao'])} · índice {s['indice_oportunidade']}</summary>",
         f"<p>{_e(s['parecer'])}</p><dl class='kv'>"]
    kv = [("PubChem", f"<a href='{_e(o.get('url_pubchem', ''))}'>CID {_e(o.get('cid'))}</a> — {_e(o.get('titulo', ''))}" if o.get("cid") else "não resolvido"),
          ("Fórmula / massa", f"{_e(o.get('formula_hill', ''))} · {_e(o.get('massa_molar', ''))} g/mol"),
          ("SMILES", _e(o.get("smiles", ""))),
          ("Átomos de N", _e("; ".join(f"{x['tipo']} [{x['basicidade']}]" for x in o.get("sitios_N", [])) or "nenhum")),
          ("Grupos ácidos", _e(_j(sorted(set(o.get("grupos_acidos", [])))) or "nenhum")),
          ("Termos buscados", _e("; ".join((r.get("termos") or {}).get("cations", []) + (r.get("termos") or {}).get("neutros", [])))),
          ("Lacunas M–X", _e(_j(s["lacunas_MX"]) or "—")),
          ("Fontes", " · ".join(f"{_e(k)}: {'ok' if v.get('ok') else '<span class=warn>falhou</span>'}" for k, v in (r.get("fontes") or {}).items()))]
    p += [f"<dt>{k}</dt><dd>{v}</dd>" for k, v in kv]
    p.append("</dl>")
    for w in o.get("avisos", []):
        p.append(f"<p class='warn small'>{_e(w)}</p>")
    if r.get("compostos"):
        p.append("<h3>Compostos reportados (amostras)</h3><div class='scroll'><table><tr><th>Composto</th><th>Categoria</th>"
                 "<th>Confiança</th><th class='num'>Depósitos</th><th class='num'>Publicações</th><th>Anos</th>"
                 "<th>Grupo espacial</th><th>Arquitetura</th></tr>")
        for g in r["compostos"]:
            yrs = f"{g['primeiro_ano']}–{g['ultimo_ano']}" if g["primeiro_ano"] else ""
            p.append(f"<tr><td>{_e(g['composto'])}</td><td>{_e(g['categoria'])}</td><td>{_e(g['melhor_confianca'])}</td>"
                     f"<td class='num'>{g['n_depositos']}</td><td class='num'>{g['n_publicacoes']}</td><td>{_e(yrs)}</td>"
                     f"<td>{_e(_j(g['grupos_espaciais']))}</td><td>{_e(_j(g['classificacoes'], ' | '))}</td></tr>")
        p.append("</table></div>")
    if r.get("estruturas"):
        p.append("<h3>Estruturas cristalinas</h3><div class='scroll'><table><tr><th>Entrada</th><th>Confiança</th>"
                 "<th>Fórmula</th><th>Arquitetura</th><th>Grupo esp.</th><th>Ano</th><th>DOI</th></tr>")
        for x in r["estruturas"][:60]:
            p.append(f"<tr><td><a href='{_e(x.get('url', ''))}'>{_e(x['fonte'])} {_e(x['id'])}</a></td>"
                     f"<td>{_e(x['confianca'])}</td><td>{_e(x.get('formula_txt', ''))}</td>"
                     f"<td>{_e(x.get('classificacao', ''))}</td><td>{_e(x.get('grupo_espacial', ''))}</td>"
                     f"<td>{_e(x.get('ano') or '')}</td><td>{_doi_link(x.get('doi'))}</td></tr>")
        p.append("</table></div>")
        if len(r["estruturas"]) > 60:
            p.append(f"<p class='small muted'>Mostrando 60 de {len(r['estruturas'])}; lista completa no Excel.</p>")
    shown = [x for x in r.get("literatura", []) if x.get("relevante")]
    if shown:
        p.append("<h3>Literatura (título/resumo citam o orgânico em contexto de haleto metálico)</h3><div class='scroll'>"
                 "<table><tr><th>Uso</th><th>Ano</th><th>Artigo</th><th class='num'>Cit.</th><th>Evidência</th></tr>")
        for x in shown[:40]:
            ev = "; ".join(x.get("pistas_estruturais", [])[:3] + x.get("pistas_aditivo", [])[:2])
            link = x.get("url") or (f"https://doi.org/{x['doi']}" if x.get("doi") else "")
            title = f"<a href='{_e(link)}'>{_e(x.get('titulo', ''))}</a>" if link else _e(x.get("titulo", ""))
            extra = f"<br><span class='small muted'>“{_e(x['trecho'][:260])}”</span>" if x.get("trecho") else ""
            p.append(f"<tr><td>{_e(x.get('uso', ''))}</td><td>{_e(x.get('ano') or '')}</td>"
                     f"<td>{title} <span class='muted'>{_e(x.get('revista', ''))}</span>{extra}</td>"
                     f"<td class='num'>{_e(x.get('citacoes') if x.get('citacoes') is not None else '')}</td>"
                     f"<td class='small'>{_e(ev)}</td></tr>")
        p.append("</table></div>")
        if len(shown) > 40:
            p.append(f"<p class='small muted'>Mostrando 40 de {len(shown)}; lista completa no Excel.</p>")
    if s["n_nao_verificaveis"]:
        p.append(f"<p class='small muted'>{s['n_nao_verificaveis']} artigo(s) retornados pelas buscas não citam o orgânico "
                 "no título/resumo (provável menção só no texto completo) — listados na aba Literatura do Excel.</p>")
    if r.get("depositos_csd"):
        nums = ", ".join(str(d["ccdc"]) for d in r["depositos_csd"] if d.get("ccdc"))
        p.append(f"<p class='small'><b>Depósitos CCDC ligados aos artigos estruturais:</b> {_e(nums)}</p>")
    p.append("</details>")
    return "".join(p)
