# Incremental CRC hypothesis-generation experiments

The experiment preserves the original notebooks, datasets, saved results and
the notebook's model architecture. It reconstructs compact inputs rather than
exporting a second full graph. CPU execution uses two threads; no API keys or
agent integration are involved.

Dependencies: pandas, NumPy, SciPy, scikit-learn and PyTorch (available through
`requirements/research.txt` plus `torch` and `scikit-learn`, or the existing
modeling environment). HGNC must be present locally at
`data/external/hgnc/hgnc_complete_set.txt` or supplied via `--hgnc`.
The reference used for this run was downloaded from
`https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt`.
The runner itself never downloads or overwrites source datasets.

```sh
python3 analysis/crc_experiments.py
python3 -m unittest discover -s tests/research -p 'test_crc_experiments.py' -v
```

Outputs use a new UTC timestamped directory under `reports/crc_experiments/`.
An explicit `--output` must be nonexistent. `--epochs` and `--seeds` allow smaller
debug runs; defaults are 40 epochs and initialization seeds42/43/44, with a
single fixed disease split seed42. Keep default settings for comparison with
the reported results. Legacy and no-PPI checks use seed42 only. Weighted subtype
profiles are inference-only interventions because this architecture has no
feedback from subtype aggregation into supervised disease training.

Read `crc_findings.md` first; `crc_experiment_report.md` gives longer methods.
The experiment **does not** establish GNN superiority or clinical efficacy.
Do not interpret deterministic ordering within exactly tied compound groups as
an experimentally meaningful difference.

Artifacts include evaluation metrics, disease split IDs, fixed test/validation
pairs, training curves, profile gene evidence, top rankings, full compressed
score arrays, pairwise rank/score diagnostics, embedding similarities, variance
decomposition, and tie-inclusive candidates. `provenance.json` fingerprints all
protected original notebooks/reports/artifacts and input datasets before and
after execution. Historical notebook outputs are not overwritten.

To regenerate the concise summary without fitting models:

```sh
python3 analysis/crc_report.py reports/crc_experiments/RUN_DIRECTORY
```

Only newly generated summary/audit files in that run directory are written.
The summary additionally consults existing compatibility drug-target and CRC
indication audits to annotate target signatures and broad CRC indication flags.
It never treats these flags as subtype response labels or evidence of literature
novelty.
