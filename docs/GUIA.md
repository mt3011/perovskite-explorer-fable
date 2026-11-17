# Guia de uso e metodologia

Este guia tem duas partes: **como usar** (seções 1–3) e **como o programa chega às respostas**
(seções 4–9), para que você possa julgar a qualidade de cada informação devolvida.

---

## 1. Usando o programa

```bash
python -m pip install -r requirements.txt
```

```bash
python -m perovskite_explorer diagnostico
```

```bash
python -m perovskite_explorer buscar minha_lista.txt -o resultados --email seu@email
```

- `diagnostico` faz uma consulta mínima a cada fonte e mostra quais estão respondendo.
- `buscar` processa a lista. O `--email` é opcional: é enviado ao Crossref/OpenAlex apenas para
  identificar o uso e dar acesso à fila mais rápida.
- **Lista**: um orgânico por linha. Aceita nome (em inglês), número CAS, `SMILES:c1ccncc1`,
  InChIKey ou `CID:1049`. Linhas com `#` são comentários.

Opções mais úteis de `buscar`:

| Opção | Efeito |
|---|---|
| `--metais Pb,Sn,Bi` | restringe os metais (padrão: Pb, Sn, Ge, Bi, Sb, In, Tl, Cu, Ag, Au, Cd, Hg, Zn, Mn, Te, Pd, Pt, Zr, Hf) |
| `--metais-extras` | acrescenta Fe, Co, Ni, Cr, V, Ti, Ga, Se, Mo, W, Re, Ru, Rh, Ir, Eu, Yb |
| `--halogenios Cl,Br,I,F` | padrão: Cl, Br, I |
| `--max-cifs 1000` | verifica mais CIFs quando há muitas candidatas só por fórmula (padrão: 300) |
| `--fontes europepmc,crossref` | escolhe as fontes de literatura (padrão: openalex, europepmc, crossref) |
| `--openalex-sem-chave` | tenta o OpenAlex com a cota anônima (pequena e compartilhada por IP; pode falhar) |
| `--offline` | usa só o que já está no cache (sem rede) |
| `--refazer` / `--atualizar` | ignora checkpoints / ignora o cache |

Outros comandos:

```bash
python -m perovskite_explorer relatorio -o resultados
```

```bash
python -m perovskite_explorer analisar-cif estrutura.cif --organico "Phenazine"
```

`relatorio` regenera os arquivos de saída sem acessar a rede. `analisar-cif` aplica a análise
estrutural (seção 6) a qualquer CIF seu — por exemplo, um baixado do site do CCDC.

**Tempo**: cerca de 3–8 minutos por orgânico na primeira vez (o COD é consultado para cada par
metal–haleto). Orgânicos com os mesmos elementos reaproveitam o cache e ficam bem mais rápidos.

**OpenAlex (recomendado)**: é a fonte de literatura com maior cobertura, mas exige uma chave de API
— gratuita, criada em https://openalex.org. Sem chave ele fica **desativado** (aparece como "não
usada" na aba Fontes) e a literatura vem só do Europe PMC e do Crossref. Para ativar, defina a
variável de ambiente antes de rodar:

```bash
export OPENALEX_API_KEY="sua_chave"
```

No Windows (PowerShell): `$env:OPENALEX_API_KEY="sua_chave"`. Também é possível gravar a chave em
um arquivo `config_local.json` (`{"openalex_api_key": "sua_chave"}`), que o git ignora.

## 2. Se a execução for interrompida

Rode **o mesmo comando** de novo. Nada do que já foi obtido é refeito:

| Proteção | Como funciona |
|---|---|
| Cache HTTP | toda resposta é gravada em `cache/http_cache.sqlite` no instante em que chega |
| Checkpoint por orgânico | cada orgânico concluído vira um JSON em `resultados/checkpoints/` (escrita atômica: nunca fica um arquivo pela metade) |
| Relatórios parciais | os arquivos de saída são regravados a cada orgânico concluído |
| Repetição automática | erros de rede e de servidor são repetidos até 5 vezes, com espera crescente |
| Isolamento de falhas | se uma fonte cai, as outras continuam; o orgânico é marcado "busca incompleta" e refeito na próxima execução |
| Disjuntor | um servidor fora do ar deixa de ser consultado naquela execução, em vez de travar por horas |

Um orgânico com **busca incompleta nunca deve ser lido como "inexplorado"**: o relatório avisa em
destaque quando isso acontece.

## 3. Lendo os resultados

| Arquivo | Conteúdo |
|---|---|
| `relatorio.html` | visão geral e ficha de cada orgânico (abra no navegador) |
| `resultados.xlsx` | todas as evidências, em abas (abaixo) |
| `resumo.csv`, `compostos.csv`, `estruturas.csv`, `literatura.csv` | as mesmas tabelas em CSV (separador `;`) |
| `dados_completos.json` | tudo, para reprocessamento |
| `execucao.log` | registro da execução |

Abas do Excel: **Resumo** (uma linha por orgânico), **Compostos** (amostras reportadas, agrupando
depósitos repetidos), **Estruturas** (uma linha por entrada cristalográfica), **Literatura** (uma
linha por artigo, com as pistas que levaram à classificação), **Depósitos CSD por artigo**,
**Descartadas** (o que foi rejeitado e por quê), **Fontes** (o que funcionou) e **Parâmetros**.

**Nível de exploração** (coluna principal do Resumo). Ele mede o uso do orgânico **como cátion na
rede de um haleto metálico**, que é o papel usual em perovskitas híbridas:

| Nível | Significado |
|---|---|
| INEXPLORADO | nenhuma estrutura e nenhum artigo em contexto de haleto metálico |
| SEM USO ESTRUTURAL CONHECIDO | só aparece como aditivo/passivante ou menção |
| USADO SÓ COMO LIGANTE OU NEUTRO | já cristalizou com haletos metálicos, mas sempre coordenado ao metal ou como molécula neutra — nunca como cátion |
| INDÍCIOS NÃO CONFIRMADOS | há evidência fraca (ver seção 5) — confira manualmente antes de considerar inédito |
| POUCO EXPLORADO | 1–2 publicações como cátion, ou um único grupo de pesquisa |
| MODERADAMENTE EXPLORADO | 3–9 publicações como cátion |
| BEM ESTABELECIDO | 10 ou mais |

O total de publicações estruturais em qualquer papel aparece em colunas separadas. Exemplo real: a
nicotinamida tem dezenas de estruturas com haletos metálicos, quase todas como ligante coordenado;
como cátion (3-carbamoilpiridínio) são poucas.

---

## 4. Metodologia: o caminho de cada orgânico

```
nome da lista
   │ 1. PubChem ........ identidade: fórmula, grafo molecular, sinônimos, sítios básicos
   │ 2. Nomes .......... variantes catiônicas (piridina → piridínio)
   ├─ 3. COD ........... busca por ELEMENTOS e por nome → fórmula → CIF → grafo + arquitetura
   ├─ 4. CSD ........... busca pelo NOME QUÍMICO do depósito (metadados públicos DataCite)
   ├─ 5. Literatura .... Europe PMC + OpenAlex + Crossref → relevância → tipo de uso
   │ 6. Cruzamento ..... artigos ↔ estruturas ↔ depósitos CCDC de cada artigo
   ▼ 7. Resumo ......... nível de exploração, o quanto foi reportado, índice de oportunidade
```

**1. Identidade (PubChem).** O nome é resolvido para um registro (CID). Do registro 2D vêm a
fórmula e o **grafo molecular** (átomos e ligações), que é a referência usada para reconhecer o
orgânico dentro dos cristais. Se a entrada for um sal, fica só o fragmento orgânico.

**2. Variantes de nome.** Num híbrido o orgânico quase sempre está protonado e muda de nome.
As variantes vêm de (a) formas catiônicas registradas no PubChem e (b) regras de nomenclatura:
`-amine → -ammonium/-aminium`, `-ine → -inium`, `-amide → -amidium`, `-ole → -olium`,
`amino- → ammonio-`, `n-aminobenzoic acid → n-carboxyanilinium`,
`pyridine-n-carboxylic acid → n-carboxypyridinium`. As regras só são aplicadas ao nome da lista,
ao nome principal/IUPAC do PubChem e a sinônimos sistemáticos, para não gerar por acidente o nome
de outro composto. Os termos usados ficam registrados no Resumo.

**3. COD.** Para cada par metal–haleto são pedidas todas as entradas que contêm os elementos do
orgânico + metal + haleto (não depende de como os autores nomearam o composto), além de buscas
textuais pelos nomes. Cada entrada passa por três testes independentes:

- **Fórmula**: `fórmula do cristal = k × (orgânico ± prótons) + metais + haletos (+ água)`, com
  balanço de H coerente (no máximo um próton por sítio básico).
- **Nome**: algum termo do orgânico no nome químico ou no título, como nome próprio
  (`dibenzo[a,c]phenazine` e `N-methylphenazinium` não contam como `phenazine`).
- **Grafo**: o CIF é baixado, a simetria é expandida, as ligações são reconstruídas por distância
  e procura-se uma molécula com **a mesma conectividade de átomos pesados** do orgânico.

Isômeros têm a mesma fórmula (fenazina e 1,10-fenantrolina são ambas C12H8N2), por isso **a
fórmula sozinha nunca conta como evidência**: ela só seleciona quais CIFs verificar.

**4. CSD sem licença.** O CCDC publica no DataCite os metadados de cada depósito do CSD, incluindo
o **nome químico sistemático** dado pelos editores do CSD, o número CCDC e o artigo associado. O
programa pesquisa esses nomes pelos termos do orgânico, exige haleto como ligante ou ânion
(`tetrabromo-`, `(mu-iodo)`, `chloride`; o cloro de `2-chloro-benzoate` não conta) e um metal da
lista. O prefixo `catena-` indica rede polimérica. **Limite**: por aqui não vêm fórmula, cela nem
coordenadas — o link do depósito leva à página do CCDC, onde a estrutura pode ser vista e o CIF
baixado gratuitamente, e então analisado com `analisar-cif`.

**5. Literatura.** A consulta é `(nomes do orgânico) E (perovskite, lead iodide, halometallate…)`.
Cada artigo retornado é reavaliado localmente e só conta se **o título ou o resumo** citarem o
orgânico como nome próprio **e** tiverem contexto de haleto metálico. O tipo de uso sai de um
placar de pistas (seção 7). Registros repetidos (mesmo DOI, ou preprint + artigo com o mesmo
título) são unidos.

**6. Cruzamento.** O artigo de cada estrutura confirmada entra na lista de literatura como
estrutural. Para os artigos estruturais, são listados todos os depósitos CCDC do artigo e marcado
se o nome químico de cada um cita o orgânico.

## 5. Níveis de confiança de uma estrutura

| Confiança | Critério | Conta como evidência? |
|---|---|---|
| **confirmada (grafo no CIF)** | molécula com a mesma conectividade encontrada no CIF | sim |
| **provável (fórmula + nome)** | fórmula compatível e nome citado; grafo inconclusivo (ex.: orgânico desordenado) | sim |
| **provável (nome sistemático no CSD)** | nome do CSD cita o orgânico; sem CIF | sim |
| possível (só fórmula) | fórmula compatível, CIF analisado mas inconclusivo | não — vira "indício" |
| possível (só nome) | nome cita, mas o grafo do CIF não confirmou | não — vira "indício" |
| descartada | isômero, ou parte orgânica ordenada sem a molécula | não (aba Descartadas) |

Ordem de verificação dos CIFs: sempre as entradas com nome citado e as de fórmula exata com o
orgânico protonado; depois as de fórmula exata com orgânico neutro, até `--max-cifs`. Se o limite
for atingido, o relatório informa quantas candidatas ficaram **pendentes**.

Cada estrutura também recebe o **papel do orgânico**: *cátion* (H a mais que a molécula neutra, ou
nome em `-ium`), *ligante coordenado ao metal* (ligação M–N/O/S medida no CIF) ou *neutro*
(co-cristal/solvato). É o papel de cátion que define o nível de exploração; quando o papel não pôde
ser determinado, a estrutura conta como cátion (lado conservador: não infla a novidade).

Só entram estruturas com **metal da lista-alvo e haleto pedido** (`--metais`, `--halogenios`);
compostos apenas com outros metais (lantanídeos, Co, Ni…) ficam de fora, a menos que você amplie a
lista.

## 6. Como a arquitetura estrutural é determinada (a partir do CIF)

1. **Ligações M–X**: metal e haleto ligados se a distância for menor que a soma dos raios
   covalentes + 0,6 Å (Pb–I até 3,45 Å). Cs, Rb, Ba e Sr são tratados como cátions do sítio A.
2. **Dimensionalidade**: percorre-se a rede M–X anotando em qual cela cada átomo é alcançado; os
   ciclos que fecham em outra cela dão vetores de translação, e o posto desses vetores é a
   dimensionalidade (0D, 1D, 2D ou 3D).
3. **Compartilhamento**: nº de haletos em ponte entre dois metais — 1 = vértices, 2 = arestas,
   3 = faces.
4. **Esfera de coordenação**: `MX6` = seis haletos; `MX4L2` = quatro haletos e dois átomos N/O/S.

"**Tipo perovskita**" só é atribuído a octaedros MX6 puros que compartilham apenas vértices
(3D, camadas 2D ou cadeias 1D); octaedros MX6 ligados por arestas/faces são "perovskitoides".
O procedimento foi conferido em estruturas de referência: CsPbBr3 (3D), Cs2PbI2Cl2 (2D),
CsNiCl3 (1D por faces) e Cs2SnI6 (0D).

## 7. Como um artigo é classificado

| Pistas de uso **estrutural** | Peso | Pistas de uso **não estrutural** | Peso |
|---|---|---|---|
| cristalografia (crystal structure, single-crystal, space group…) | 3 | aditivo/passivação (passivat-, additive, defect, interface, dopant…) | 3 |
| orgânico na forma catiônica (`-ium`) | 3 | dispositivo (solar cell, PCE, fill factor…) | 2 |
| orgânico definido em fórmula ("PZ = phenazine") | 3 | processamento/camadas (HTL, ETL, antisolvent, SAM…) | 2 |
| dimensionalidade (2D, layered, Ruddlesden…) | 2 | | |
| arquitetura (octahedra, chains, spacer, templat-…) | 2 | | |
| halometalato / híbrido orgânico–inorgânico | 2 | | |
| química de coordenação (complex, ligand…) | 2 | | |
| propriedades do material; síntese de novo composto; "structur-" | 1 | | |

Pistas fora da frase que cita o orgânico valem metade. **Estrutural**: ≥ 4 pontos e pelo menos 2
a mais que o outro lado (o mesmo vale, invertido, para **aditivo**). Empate com ≥ 4: **ambíguo**.
Abaixo disso: **menção**. Um artigo com estrutura cristalina confirmada é sempre estrutural.
Todas as pistas e o trecho do resumo ficam na aba Literatura para conferência.

## 8. "O quanto já foi reportado" e o índice de oportunidade

- **Amostras únicas**: estruturas do mesmo artigo com os mesmos metais e haletos costumam estar
  nos dois bancos; para cada (artigo, M, X) conta-se o maior dos dois totais (COD ou CSD).
- **Publicações como cátion**: DOIs distintos de estruturas em que o orgânico é cátion, mais os
  artigos classificados como estruturais cuja estrutura não foi localizada nos bancos.
- **Publicações estruturais (total)**: o mesmo, incluindo ligante coordenado e forma neutra.
- **Grupos independentes** (aproximado): artigos ligados por coautoria (2 ou mais autores em
  comum) contam como um grupo. Um composto com 5 artigos do mesmo grupo foi pouco reproduzido.
- **Combinações M–X** e **lacunas**: quais pares metal–haleto já têm estrutura com o orgânico e
  quais do núcleo (Pb, Sn, Bi, Sb, Cu, Ag × Cl, Br, I) não foram encontrados.

**Índice de oportunidade = novidade × viabilidade** (0–100; serve só para ordenar a lista):

| Novidade | | Viabilidade (formar cátion) | |
|---|---|---|---|
| INEXPLORADO | 100 | amina alifática ou carga permanente | 1,0 |
| SEM USO ESTRUTURAL | 90 | N sp2 básico (piridina, imidazol, amidina) | 0,85 |
| SÓ LIGANTE/NEUTRO | 80 | só base fraca (tipo anilina) | 0,7 |
| INDÍCIOS | 75 | nenhum sítio básico | 0,3 |
| POUCO | 65 | não resolvido no PubChem | 0,5 |
| MODERADO / ESTABELECIDO | 35 / 10 | | |

A classificação dos nitrogênios é por grupo funcional, lida do grafo: amida, sulfonamida, nitrila,
nitro e N tipo pirrol **não** contam como básicos.

## 9. Limitações (leia antes de concluir que algo é inédito)

- **Nenhuma busca prova ausência.** "INEXPLORADO" significa "não encontrado nestas fontes com
  estes nomes". Confira os *termos de busca* no Resumo: se o cátion tiver um nome usual que não
  está lá, acrescente-o à lista como uma linha extra.
- **CSD por nome** depende do nome sistemático: abreviações não existem no CSD, mas um nome
  sistemático diferente do esperado pode escapar. A busca por subestrutura (seção abaixo)
  elimina essa limitação.
- **COD** é bem menor que o CSD (que tem mais de 1,2 milhão de depósitos); muitas estruturas só existem no CSD.
- **Literatura**: só título e resumo são lidos. Artigos em que o orgânico aparece apenas no texto
  completo, ou sob abreviação, ficam como "não verificáveis" (listados, mas sem contar).
  O classificador por pistas erra em casos de fronteira — por isso as pistas ficam visíveis.
- **Desordem**: se a parte orgânica do CIF está desordenada, o grafo pode não fechar; a entrada
  fica como "provável" (se o nome confirmar) ou "possível".
- **Arquitetura** usa limites de distância fixos: ligações secundárias longas (Sb, Bi, Sn com par
  isolado ativo) podem ficar de fora e reduzir a dimensionalidade aparente.
- **Escopo**: haletos (Cl, Br, I; F opcional). Pseudo-haletos (SCN, formiato, azida) e
  calcogenetos não são procurados.

## Quando você tiver a API do CSD

O conector `perovskite_explorer/csd.py` já está escrito para a **CSD Python API** (`ccdc`): faz
busca por **subestrutura** com o esqueleto do orgânico, exige um halogênio e entrega cada acerto,
com o CIF, à mesma verificação de grafo e análise estrutural usadas no COD. Ele é ativado
sozinho quando o pacote `ccdc` estiver instalado (o `diagnostico` mostra o estado).
**Ele ainda não foi testado**, por falta de licença: foi escrito a partir da documentação pública.
Na primeira execução com a licença, compare o resultado de um orgânico conhecido (ex.: fenazina)
com o que a busca por nome já devolve. Se o seu acesso for por outra interface, basta reescrever
a função `search` desse arquivo mantendo o formato de retorno.
