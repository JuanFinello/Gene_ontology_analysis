# Gene Ontology enrichment analysis

GO enrichment of microarray data from the moss *Physcomitrium patens*, comparing a wild
type line against a *pten* knockout. Research assistant work at the Plant Physiology
Laboratory, National University of Córdoba.

## Data

A microarray table of 32,559 genes measured in two conditions, wild type and *pten*
knockout, with three replicates each.

## Pipeline

| Script | What it does |
|---|---|
| `1_preparacion_de_datos.R` | Per-gene t-test between the two conditions, producing a fold change and a p-value for every gene |
| `2_construcción_del_objeto_GO.R` | Builds the `topGOdata` objects: gene universe, selection function and GO annotation mapping |
| `3_analisis_de_enriquecimiento.R` | Fisher exact test over the GO graph, p-value distribution and a ranked table of enriched terms |

Enrichment is run separately per GO ontology (cellular component, biological process,
molecular function) and per direction of change.

## Stack

R · Bioconductor · topGO
