# Perovskite Explorer

Recebe uma lista de compostos orgânicos (um `.txt`) e responde, para cada um:

- **Já foi usado estruturalmente** em haletos metálicos híbridos (perovskitas 3D/2D, cadeias 1D,
  unidades 0D, perovskitoides, compostos de coordenação com haleto)?
- **Qual é a estrutura** das amostras reportadas (dimensionalidade, tipo de poliedro, grupo espacial,
  cela, papel do orgânico)?
- **O quanto foi reportado** (nº de estruturas, publicações, grupos independentes, citações, anos)?
- **Vale a pena estudar?** Nível de exploração + índice de oportunidade, com todas as evidências
  rastreáveis (links para o banco e para o artigo).

Funciona **sem nenhuma IA/LLM**: só consultas a bancos públicos e regras determinísticas.

| Fonte | Para quê |
|---|---|
| PubChem | identidade do orgânico, fórmula, grafo molecular, sinônimos, sais conhecidos |
| COD (Crystallography Open Database) | estruturas com CIF: confirmação por grafo e análise estrutural |
| CSD/CCDC via DataCite | busca no CSD pelo nome químico e depósitos de cada artigo (sem licença) |
| CSD Python API | opcional, ativado automaticamente quando você tiver a licença |
| Europe PMC, Crossref, OpenAlex | literatura (título/resumo), citações |

## Instalação

```bash
git clone https://github.com/mt3011/perovskite-explorer-fable.git
cd perovskite-explorer-fable
python -m pip install -r requirements.txt
```

Requer Python 3.9 ou mais novo e acesso à internet.

## Uso rápido

```bash
python -m perovskite_explorer diagnostico
```

```bash
python -m perovskite_explorer buscar exemplos/organicos.txt -o resultados
```

Abra `resultados/relatorio.html` (visão geral) e `resultados/resultados.xlsx` (todas as evidências).
Se a execução for interrompida, **rode o mesmo comando de novo**: ela continua de onde parou.

A lista é um arquivo de texto com um orgânico por linha (nome em inglês, CAS, `SMILES:...`,
InChIKey ou `CID:1234`); linhas com `#` são comentários. Veja `exemplos/organicos.txt`.

## Documentação

**[docs/GUIA.md](docs/GUIA.md)** — como usar, como ler os resultados e toda a metodologia de busca
(critérios, limites e o que foi validado).

## Testes

```bash
python -m pip install pytest
python -m pytest
```
