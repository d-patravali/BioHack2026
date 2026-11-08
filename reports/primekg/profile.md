# Supplied PrimeKG dataset profile

Raw rows: 8,100,498. Unique nodes: 129,375.

## Node types

| Type | Nodes |
| --- | ---: |
| anatomy | 14,035 |
| biological_process | 28,642 |
| cellular_component | 4,176 |
| disease | 17,080 |
| drug | 7,957 |
| effect/phenotype | 15,311 |
| exposure | 818 |
| gene/protein | 27,671 |
| molecular_function | 11,169 |
| pathway | 2,516 |

## Typed relationships (raw row counts)

| Relation | Display label | X type | Y type | Rows |
| --- | --- | --- | --- | ---: |
| anatomy_anatomy | parent-child | anatomy | anatomy | 28,064 |
| anatomy_protein_absent | expression absent | anatomy | gene/protein | 19,887 |
| anatomy_protein_absent | expression absent | gene/protein | anatomy | 19,887 |
| anatomy_protein_present | expression present | anatomy | gene/protein | 1,518,203 |
| anatomy_protein_present | expression present | gene/protein | anatomy | 1,518,203 |
| bioprocess_bioprocess | parent-child | biological_process | biological_process | 105,772 |
| bioprocess_protein | interacts with | biological_process | gene/protein | 144,805 |
| bioprocess_protein | interacts with | gene/protein | biological_process | 144,805 |
| cellcomp_cellcomp | parent-child | cellular_component | cellular_component | 9,690 |
| cellcomp_protein | interacts with | cellular_component | gene/protein | 83,402 |
| cellcomp_protein | interacts with | gene/protein | cellular_component | 83,402 |
| contraindication | contraindication | disease | drug | 30,675 |
| contraindication | contraindication | drug | disease | 30,675 |
| disease_disease | parent-child | disease | disease | 64,388 |
| disease_phenotype_negative | phenotype absent | disease | effect/phenotype | 1,193 |
| disease_phenotype_negative | phenotype absent | effect/phenotype | disease | 1,193 |
| disease_phenotype_positive | phenotype present | disease | effect/phenotype | 150,317 |
| disease_phenotype_positive | phenotype present | effect/phenotype | disease | 150,317 |
| disease_protein | associated with | disease | gene/protein | 80,411 |
| disease_protein | associated with | gene/protein | disease | 80,411 |
| drug_drug | synergistic interaction | drug | drug | 2,672,628 |
| drug_effect | side effect | drug | effect/phenotype | 64,784 |
| drug_effect | side effect | effect/phenotype | drug | 64,784 |
| drug_protein | carrier | drug | gene/protein | 864 |
| drug_protein | carrier | gene/protein | drug | 864 |
| drug_protein | enzyme | drug | gene/protein | 5,317 |
| drug_protein | enzyme | gene/protein | drug | 5,317 |
| drug_protein | target | drug | gene/protein | 16,380 |
| drug_protein | target | gene/protein | drug | 16,380 |
| drug_protein | transporter | drug | gene/protein | 3,092 |
| drug_protein | transporter | gene/protein | drug | 3,092 |
| exposure_bioprocess | interacts with | biological_process | exposure | 1,625 |
| exposure_bioprocess | interacts with | exposure | biological_process | 1,625 |
| exposure_cellcomp | interacts with | cellular_component | exposure | 10 |
| exposure_cellcomp | interacts with | exposure | cellular_component | 10 |
| exposure_disease | linked to | disease | exposure | 2,304 |
| exposure_disease | linked to | exposure | disease | 2,304 |
| exposure_exposure | parent-child | exposure | exposure | 4,140 |
| exposure_molfunc | interacts with | exposure | molecular_function | 45 |
| exposure_molfunc | interacts with | molecular_function | exposure | 45 |
| exposure_protein | interacts with | exposure | gene/protein | 1,212 |
| exposure_protein | interacts with | gene/protein | exposure | 1,212 |
| indication | indication | disease | drug | 9,388 |
| indication | indication | drug | disease | 9,388 |
| molfunc_molfunc | parent-child | molecular_function | molecular_function | 27,148 |
| molfunc_protein | interacts with | gene/protein | molecular_function | 69,530 |
| molfunc_protein | interacts with | molecular_function | gene/protein | 69,530 |
| off-label use | off-label use | disease | drug | 2,568 |
| off-label use | off-label use | drug | disease | 2,568 |
| pathway_pathway | parent-child | pathway | pathway | 5,070 |
| pathway_protein | interacts with | gene/protein | pathway | 42,646 |
| pathway_protein | interacts with | pathway | gene/protein | 42,646 |
| phenotype_phenotype | parent-child | effect/phenotype | effect/phenotype | 37,472 |
| phenotype_protein | associated with | effect/phenotype | gene/protein | 3,330 |
| phenotype_protein | associated with | gene/protein | effect/phenotype | 3,330 |
| protein_protein | ppi | gene/protein | gene/protein | 642,150 |

## Phenotype–disease–gene extraction

60,062 nodes; 1,214,512 raw rows; 1,214,512 unique ordered edges.
Removed 0 repeated ordered rows; multiplicity is retained in `raw_row_count`.
1,214,512 exported ordered edges have a reverse edge with identical relation labels.

See ../../docs/research/primekg_exploration.md for semantics, limitations, and construction examples.
