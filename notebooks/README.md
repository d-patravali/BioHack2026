# Graph integration and GNN experiments

Install optional dependencies in a modeling environment:

```sh
python -m pip install -r requirements/notebooks.txt
jupyter lab
```

Launch Jupyter from the repository root or a directory inside it. Each notebook
finds the checkout root before resolving data paths. Original outputs and
metadata are preserved; the reorganization did not run any cells. Review the
code before running: cells can download references, install a package, train
models, and overwrite generated files.

| Notebook | Inputs | Generated outputs |
| --- | --- | --- |
| `integration/PrimeKG_ECoGraph_merge.ipynb` | `data/raw/primekg/kg.csv`, `ECoGraph/data/coad_mut_data.csv`, `ECoGraph/data/neptune_csv/`, HGNC reference | `data/processed/primekg_ecograph/` graph tables and evidence sidecars |
| `models/colorectal/PrimeKG_ECoGraph_GNN.ipynb` | Exported graph tables | `artifacts/predictions/colorectal/subtype_drug_rankings.csv` |
| `models/breast/PrimeKG_HER2_Breast_Cancer_GNN.ipynb` | Exported graph tables; GDSC workbook | `artifacts/predictions/breast/` mapping audit, held-out predictions, and rankings |

The merge notebook downloads HGNC into `data/external/hgnc/` if absent. The
breast notebook downloads GDSC into `data/external/gdsc/` if absent. Other
existing download logic remains in the notebooks; the notebooks are not
guaranteed offline or read-only. No download was performed during relocation.

The source ECoGraph `edges_has_gene_part001.csv` is blank; the merge notebook
uses the local binary mutation matrix instead. The compatibility report explains
that distinction and the scientific limitations. Existing notebook outputs are
historical evidence of teammate runs, not a verification of the relocated code.

Review modeling assumptions—including treatment of `drug_protein` display
relations, subtype representation, and evaluation—before exposing predictions
to the agent. The agent is currently connected only to its static subtype tool.
