# CRCSC–PrimeKG subtype–drug prototype

This separate workflow exports a small typed CMS overlay on the read-only full
PrimeKG, builds a compact molecular message-passing view, and generates
**unvalidated** CMS drug-association scores. It preserves the existing agent,
notebooks, raw data and historical outputs. No literature review is automatic.

```sh
python3 -m unittest discover -s tests/research -p 'test_crcsc_drug_gnn.py' -v
python3 analysis/crcsc_drug_gnn.py
python3 analysis/crcsc_drug_report.py reports/crcsc_drug_gnn/RUN_ID
```

Use a research environment with pandas, NumPy, SciPy, scikit-learn, PyTorch and
threadpoolctl; CPU uses two threads. No external graph package or API is needed.
Default source: `reports/crcsc/20261108_172638/selected_gene_signatures.csv` plus
its completed manifest. Source genes/weights were selected using discovery only.
Other inputs: raw PrimeKG and its existing node table. No ECoGraph files are used.

Defaults: 40 epochs, seeds42/43/44, matched disease split42. `--epochs`, `--seeds`,
`--source-run` and `--output` support reproducibility; output must not exist.
Review the completed run's `findings.md` and `manifest.json` before interpreting
candidate CSVs. The full base is referenced, not duplicated; overlay nodes and
edges are separately exported. Original index aliases are audited.

Models: existing PPI/target pooling, added gene–pathway–gene relation messages,
and a retrained disease-context ablation. The training target is general disease
indication recovery against sampled unknown pairs. There are no CMS therapeutic
labels. Signed expression pooling at inference is an explicit domain-transfer
assumption, not a trained drug-action-direction model. Magnitude/up/down variants
and profile nulls are diagnostics, not independently validated predictors.

Cutoff ties are preserved; identical drug-target neighborhoods remain
indistinguishable without additional compound information. The code does not
register these results as agent tools or call them clinically validated.
`literature_review_queue.csv` has explicit NOT REVIEWED evidence fields. A paper
can support/contradict an association but cannot alone verify therapeutic efficacy.
Independent subtype-sensitive experimental evidence remains necessary.
