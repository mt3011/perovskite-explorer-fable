"""Análise estrutural determinística de um CIF.

Responde a duas perguntas, sem depender de nomes:
  1. O orgânico procurado está mesmo no cristal? (isomorfismo do grafo de átomos pesados)
  2. Como é a parte inorgânica? (dimensionalidade 0D/1D/2D/3D da rede M-X, tipo de
     compartilhamento dos poliedros, número de coordenação, orgânico coordenado ao metal?)
"""
from __future__ import annotations

import itertools
from collections import Counter

import numpy as np

from .chem import (A_SITE_METALS, ALKALI_LIGHT, COVALENT_RADII, HALOGENS, METALS,
                   ORGANIC_ELEMENTS, hill_formula)

TOL_COVALENT = 0.45   # Å somados aos raios covalentes para ligações da parte orgânica
TOL_METAL_X = 0.60    # idem para M-haleto (inclui ligações longas de Pb/Sn/Bi/Sb)
TOL_METAL_L = 0.50    # idem para M-N/O/S/P (coordenação do orgânico)
MAX_ATOMS = 3000


class CifError(Exception):
    pass


def _read_block(cif_text: str):
    import gemmi

    doc = gemmi.cif.read_string(cif_text)
    for block in doc:
        if block.find_loop("_atom_site_fract_x"):
            return block
    raise CifError("CIF sem coordenadas atômicas")


def _symops(block, small):
    import gemmi

    for tag in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz",
                "_space_group_symop.operation_xyz"):
        col = block.find_loop(tag)
        if col:
            ops = []
            for raw in col:
                s = gemmi.cif.as_string(raw).replace(" ", "")
                try:
                    ops.append(gemmi.Op(s))
                except Exception:  # noqa: BLE001
                    continue
            if ops:
                return ops
    sg = small.spacegroup
    if sg is None and small.spacegroup_hm:
        sg = gemmi.find_spacegroup_by_name(small.spacegroup_hm)
    if sg is not None:
        return list(sg.operations())
    return [gemmi.Op("x,y,z")]


def _orth_matrix(cell) -> np.ndarray:
    """Matriz cujas linhas são os vetores a, b, c em coordenadas cartesianas (Å)."""
    a, b, c = cell.a, cell.b, cell.c
    al, be, ga = np.radians([cell.alpha, cell.beta, cell.gamma])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz = np.sqrt(max(c * c - cx * cx - cy * cy, 1e-12))
    return np.array([[a, 0.0, 0.0], [b * np.cos(ga), b * np.sin(ga), 0.0], [cx, cy, cz]])


def _load_atoms(cif_text: str):
    """Lê o CIF, escolhe a componente majoritária da desordem e expande a simetria."""
    import gemmi

    block = _read_block(cif_text)
    small = gemmi.make_small_structure_from_block(block)
    cell = small.cell
    if not (cell.a > 1 and cell.b > 1 and cell.c > 1):
        raise CifError("cela unitária ausente ou inválida")
    M = _orth_matrix(cell)

    raw = []
    for s in small.sites:
        el = s.element.name
        if el == "D":
            el = "H"
        if el not in COVALENT_RADII:
            continue
        occ = s.occ if s.occ > 0 else 1.0
        grp = int(getattr(s, "disorder_group", 0) or 0)
        raw.append({"el": el, "f": np.array([s.fract.x, s.fract.y, s.fract.z]), "occ": occ,
                    "grp": grp, "dis": grp != 0 or occ < 0.99})
    if not raw:
        raise CifError("nenhum átomo reconhecido")

    disorder = any(a["grp"] != 0 or (a["occ"] < 0.99 and a["el"] != "H") for a in raw)
    groups: dict = {}
    for a in raw:
        if a["grp"] != 0:
            groups.setdefault(a["grp"], []).append(a["occ"])
    if groups:
        major = max(groups, key=lambda g: (sum(groups[g]) / len(groups[g]), g > 0, -abs(g)))
        raw = [a for a in raw if a["grp"] in (0, major)]

    ops = _symops(block, small)
    den = float(gemmi.Op.DEN)
    rots = [np.array(op.rot, dtype=float) / den for op in ops]
    trans = [np.array(op.tran, dtype=float) / den for op in ops]

    # ordena por ocupação: em sobreposições fica o átomo mais ocupado
    raw.sort(key=lambda a: -a["occ"])
    els, fracs, dis = [], [], []
    acc = np.zeros((0, 3))
    for a in raw:
        for R, t in zip(rots, trans):
            f = (R @ a["f"] + t) % 1.0
            if len(acc):
                d = acc - f
                d -= np.round(d)
                dist = np.linalg.norm(d @ M, axis=1)
                if dist.min() < (0.4 if a["el"] == "H" else 0.55):
                    continue
            acc = np.vstack([acc, f])
            els.append(a["el"])
            fracs.append(f)
            dis.append(a["dis"])
            if len(els) > MAX_ATOMS:
                raise CifError(f"estrutura grande demais (> {MAX_ATOMS} átomos na cela)")

    def tag(*names):
        for n in names:
            v = block.find_value(n)
            if v and v not in ("?", "."):
                return gemmi.cif.as_string(v).strip()
        return ""

    meta = {
        "formula_cif": tag("_chemical_formula_sum"),
        "nome_cif": tag("_chemical_name_systematic", "_chemical_name_common"),
        "grupo_espacial": tag("_space_group_name_H-M_alt", "_symmetry_space_group_name_H-M")
        or (small.spacegroup_hm or ""),
        "temperatura_K": tag("_diffrn_ambient_temperature", "_cell_measurement_temperature"),
        "R": tag("_refine_ls_R_factor_gt", "_refine_ls_R_factor_obs", "_refine_ls_R_factor_all"),
        "cela": [round(x, 4) for x in (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)],
        "volume": round(cell.volume, 2),
    }
    return els, np.array(fracs), M, disorder, dis, meta


def _bond_tables(els):
    """Matriz de distâncias máximas de ligação para cada par de átomos, por classe."""
    n = len(els)
    r = np.array([COVALENT_RADII[e] for e in els])
    is_h = np.array([e == "H" for e in els])
    is_x = np.array([e in HALOGENS for e in els])
    is_m = np.array([e in METALS for e in els])
    is_org = np.array([e in ORGANIC_ELEMENTS for e in els])
    is_lig = np.array([e in ("N", "O", "S", "P") for e in els])
    light = np.array([e in ALKALI_LIGHT for e in els])
    base = r[:, None] + r[None, :]
    cut = np.zeros((n, n))
    nonmetal = is_org | is_x
    # covalentes entre não-metais (orgânico e C-X)
    cov = nonmetal[:, None] & nonmetal[None, :]
    cut[cov] = (base + TOL_COVALENT)[cov]
    cut[is_h[:, None] & is_h[None, :]] = 0.0
    cut[is_x[:, None] & is_x[None, :]] = 0.0      # X...X (poli-haletos) não entram
    cut[(is_x[:, None] & is_h[None, :]) | (is_h[:, None] & is_x[None, :])] = 0.0  # H...X é lig. de H
    mx = (is_m[:, None] & is_x[None, :]) | (is_x[:, None] & is_m[None, :])
    tol_mx = np.where(light[:, None] | light[None, :], 0.30, TOL_METAL_X)
    cut[mx] = (base + tol_mx)[mx]
    ml = (is_m[:, None] & is_lig[None, :]) | (is_lig[:, None] & is_m[None, :])
    cut[ml] = (base + TOL_METAL_L)[ml]
    return cut, is_m, is_x


def _find_bonds(els, fracs, M):
    """Todas as ligações, considerando imagens periódicas. Retorna [(i, j, (n1,n2,n3))], i<=j."""
    cut, _, _ = _bond_tables(els)
    cart = fracs @ M
    n = len(els)
    bonds = []
    for shift in itertools.product((-1, 0, 1), repeat=3):
        sv = np.array(shift, dtype=float) @ M
        for start in range(0, n, 400):
            block = cart[start:start + 400]
            d = np.linalg.norm(block[:, None, :] - (cart[None, :, :] + sv), axis=2)
            ii, jj = np.nonzero((d < cut[start:start + 400]) & (d > 0.45))
            for i, j, dist in zip(ii + start, jj, d[ii, jj]):
                # átomo j deslocado de `shift` está ligado a i; guarda cada par uma vez
                if i < j or (i == j and shift > (0, 0, 0)):
                    bonds.append((int(i), int(j), shift, float(dist)))
    return bonds


def _components(nodes, edges):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        parent[find(a)] = find(b)
    comps: dict = {}
    for n in nodes:
        comps.setdefault(find(n), []).append(n)
    return list(comps.values())


def _dimensionality(nodes, edges) -> list:
    """Dimensionalidade periódica de cada componente da rede (0 a 3).

    edges: [(u, v, shift)] onde v+shift está ligado a u. Percorre o grafo atribuindo a cada
    nó a cela em que foi alcançado; cada ciclo que fecha em outra cela dá um vetor de
    translação. O posto desses vetores é a dimensionalidade.
    """
    adj: dict = {n: [] for n in nodes}
    for u, v, s in edges:
        s = np.array(s)
        adj[u].append((v, s))
        adj[v].append((u, -s))
    seen: dict = {}
    dims = []
    for root in nodes:
        if root in seen:
            continue
        seen[root] = np.zeros(3, dtype=int)
        stack, vecs, members = [root], [], [root]
        while stack:
            u = stack.pop()
            for v, s in adj[u]:
                if v not in seen:
                    seen[v] = seen[u] + s
                    stack.append(v)
                    members.append(v)
                else:
                    cyc = seen[u] + s - seen[v]
                    if cyc.any():
                        vecs.append(cyc)
        rank = int(np.linalg.matrix_rank(np.array(vecs))) if vecs else 0
        dims.append((rank, members))
    return dims


def analyze_cif(cif_text: str, query: dict | None = None) -> dict:
    """Analisa um CIF. `query` = {"heavy_elements": [...], "heavy_edges": [[i, j], ...]}.

    Nunca lança exceção: em caso de problema devolve {"ok": False, "erro": ...}.
    """
    try:
        return _analyze(cif_text, query)
    except CifError as exc:
        return {"ok": False, "erro": str(exc)}
    except Exception as exc:  # noqa: BLE001 - CIFs reais têm todo tipo de irregularidade
        return {"ok": False, "erro": f"{type(exc).__name__}: {exc}"}


def _analyze(cif_text, query):
    import networkx as nx

    els, fracs, M, disorder, dis, meta = _load_atoms(cif_text)
    n = len(els)
    bonds = _find_bonds(els, fracs, M)
    nbrs: dict = {i: [] for i in range(n)}
    for i, j, s, d in bonds:
        nbrs[i].append((j, s, d))
        nbrs[j].append((i, tuple(-x for x in s), d))

    framework_metal = [e in METALS and e not in A_SITE_METALS for e in els]
    # haleto "inorgânico": halogênio que não está ligado covalentemente a carbono
    inorganic_x = [e in HALOGENS and not any(els[j] == "C" for j, _, _ in nbrs[i])
                   for i, e in enumerate(els)]
    org_nodes = [i for i, e in enumerate(els)
                 if (e in ORGANIC_ELEMENTS or (e in HALOGENS and not inorganic_x[i]))]
    org_set = set(org_nodes)
    org_edges = [(i, j) for i, j, s, d in bonds if i in org_set and j in org_set and i != j]
    comps = _components(org_nodes, org_edges)
    comp_of = {i: c for c, members in enumerate(comps) for i in members}

    # ---------- componentes moleculares ----------
    summary = Counter()
    organic_graphs = []
    for members in comps:
        heavy_m = [i for i in members if els[i] != "H"]
        if not heavy_m:
            continue
        f = Counter(els[i] for i in members)
        hf = Counter(els[i] for i in heavy_m)
        label = hill_formula(dict(f))
        if "C" in hf:
            g = nx.Graph()
            for i in heavy_m:
                g.add_node(i, el=els[i])
            ms = set(heavy_m)
            g.add_edges_from((i, j) for i, j in org_edges if i in ms and j in ms)
            organic_graphs.append((g, dict(hf), f.get("H", 0), members))
        elif set(hf) == {"O"} and len(heavy_m) == 1:
            label = "H2O/OH/O"
        summary[label] += 1

    # ---------- o orgânico procurado está presente? ----------
    organic = {"status": "nao_avaliado"}
    matched_members: set = set()
    if query and query.get("heavy_elements"):
        q = nx.Graph()
        for idx, e in enumerate(query["heavy_elements"]):
            q.add_node(idx, el=e)
        q.add_edges_from(tuple(e) for e in query["heavy_edges"])
        qf = dict(Counter(query["heavy_elements"]))
        same_formula, hits, h_found = 0, 0, []
        for g, hf, nh, members in organic_graphs:
            if hf != qf:
                continue
            same_formula += 1
            if nx.is_isomorphic(g, q, node_match=lambda a, b: a["el"] == b["el"]):
                hits += 1
                h_found.append(nh)
                matched_members.update(members)
        if hits:
            organic = {"status": "confirmado", "n_moleculas_na_cela": hits,
                       "h_por_molecula_no_cif": sorted(set(h_found))}
        elif not organic_graphs:
            organic = {"status": "inconclusivo", "motivo": "CIF sem átomos de carbono localizados"}
        elif same_formula:
            organic = {"status": "isomero",
                       "motivo": "há molécula com a mesma fórmula, mas conectividade diferente (isômero)"}
        elif any(dis[i] for _, _, _, members in organic_graphs for i in members if els[i] != "H"):
            # só a desordem na PARTE ORGÂNICA impede a conclusão (desordem em ânions/água não)
            organic = {"status": "inconclusivo",
                       "motivo": "parte orgânica desordenada; o grafo não pôde ser reconstruído"}
        else:
            organic = {"status": "ausente",
                       "motivo": "parte orgânica ordenada sem nenhuma molécula com a composição do orgânico"}

    # ---------- parte inorgânica ----------
    metal_idx = [i for i in range(n) if framework_metal[i]]
    mx_edges, coord_x, coord_l, bond_lengths = [], Counter(), Counter(), []
    organic_bound = False
    target_bound = False
    for i, j, s, d in bonds:
        for u, v, sh in ((i, j, s), (j, i, tuple(-x for x in s))):
            if framework_metal[u] and inorganic_x[v]:
                mx_edges.append((u, v, sh))
                coord_x[u] += 1
                bond_lengths.append(d)
            elif framework_metal[u] and els[v] in ("N", "O", "S", "P"):
                coord_l[u] += 1
                c = comp_of.get(v)
                if c is not None and any(els[k] == "C" for k in comps[c]):
                    organic_bound = True
                    if v in matched_members:
                        target_bound = True
    net_nodes = sorted({u for u, _, _ in mx_edges} | {v for _, v, _ in mx_edges})
    dims = _dimensionality(net_nodes, mx_edges) if mx_edges else []
    dim_max = max((r for r, _ in dims), default=None)

    # compartilhamento de poliedros: nº de haletos em ponte entre cada par de metais
    x_to_m: dict = {}
    for u, v, sh in mx_edges:
        # metal u vê haleto v em +sh  => visto do haleto, o metal está em -sh
        x_to_m.setdefault(v, set()).add((u, tuple(-x for x in sh)))
    pair_shared: Counter = Counter()
    for v, ms in x_to_m.items():
        for (m1, s1), (m2, s2) in itertools.combinations(sorted(ms), 2):
            rel = tuple(b - a for a, b in zip(s1, s2))
            # (m1, m2, rel) e (m2, m1, -rel) são o mesmo par; sorted() garante m1 <= m2
            key = (m1, m2, rel) if m1 < m2 else (m1, m2, max(rel, tuple(-x for x in rel)))
            pair_shared[key] += 1
    sharing = Counter()
    for cnt in pair_shared.values():
        sharing["vértices" if cnt == 1 else "arestas" if cnt == 2 else "faces"] += 1

    cn_values = sorted({coord_x[m] for m in metal_idx if coord_x[m] > 0})
    spheres = sorted({(f"MX{coord_x[m]}" if coord_x[m] else "M") + (f"L{coord_l[m]}" if coord_l[m] else "")
                      for m in metal_idx if coord_x[m] or coord_l[m]})
    metals_present = sorted({els[i] for i in metal_idx})
    a_site = sorted({e for e in els if e in A_SITE_METALS})
    halides = sorted({els[i] for i in range(n) if inorganic_x[i]})

    result = {
        "ok": True,
        "n_atomos_cela": n,
        "desordem": bool(disorder),
        "organico": organic,
        "componentes": dict(summary),
        "metais": metals_present,
        "cations_inorganicos_sitio_A": a_site,
        "haletos": halides,
        "dimensionalidade": dim_max,
        "compartilhamento": dict(sharing),
        "numeros_de_coordenacao_MX": cn_values,
        "esferas_de_coordenacao": spheres,
        "distancia_MX_media": round(float(np.mean(bond_lengths)), 3) if bond_lengths else None,
        "organico_coordenado_ao_metal": organic_bound,
        "organico_alvo_coordenado_ao_metal": target_bound,
        "ligantes_NOS_por_metal_max": max(coord_l.values(), default=0),
        **meta,
    }
    result["classificacao"] = classify(result)
    return result


def classify(r: dict) -> str:
    """Rótulo legível da arquitetura estrutural.

    Notação da esfera de coordenação: MX6 = metal ligado a 6 haletos; MX4L2 = 4 haletos e
    2 átomos N/O/S/P (L). "Tipo perovskita" só é atribuído a octaedros MX6 puros que
    compartilham apenas vértices.
    """
    if not r.get("metais"):
        if r.get("haletos"):
            return "Sal orgânico de haleto, sem metal (precursor ou composto livre de metal)"
        return "Estrutura sem metal e sem haleto iônico"
    spheres = r.get("esferas_de_coordenacao") or []
    mx = f"M = {', '.join(r['metais'])}; X = {', '.join(r.get('haletos') or []) or '—'}"
    pure = bool(spheres) and all("L" not in sp for sp in spheres)
    units = "/".join("octaedros MX6" if sp == "MX6" else f"unidades {sp}" for sp in spheres) or "unidades M–X"
    dim = r.get("dimensionalidade")
    sh = r.get("compartilhamento") or {}
    if dim is None:
        base = "Sem ligação M–X: haleto fora da esfera do metal (contra-íon)"
        if spheres:
            base += f"; esfera {'/'.join(spheres)}"
    elif dim == 0:
        base = (f"0D — clusters isolados de {units} (compartilhando {_join(sh)})" if sh
                else f"0D — {units} isolados(as)")
    else:
        name = {1: "1D — cadeias", 2: "2D — camadas", 3: "3D — rede tridimensional"}[dim]
        base = f"{name} de {units} compartilhando {_join(sh) or '?'}"
        if pure and spheres == ["MX6"]:
            base += " (tipo perovskita)" if set(sh) == {"vértices"} else " (perovskitoide)"
    if r.get("organico_alvo_coordenado_ao_metal"):
        base += "; o orgânico procurado está COORDENADO ao metal (composto de coordenação)"
    elif r.get("organico_coordenado_ao_metal"):
        base += "; outro orgânico coordenado ao metal"
    elif not pure and spheres:
        base += "; metal também ligado a O/N de água/solvente/oxo"
    return f"{base} [{mx}]"


def _join(sh: dict) -> str:
    order = [k for k in ("vértices", "arestas", "faces") if k in sh]
    return " e ".join(order)
