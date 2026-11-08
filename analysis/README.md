# ECoGraph–PrimeKG compatibility

For the separate CRC GNN sensitivity/ablation experiment, see
[CRC_EXPERIMENTS.md](CRC_EXPERIMENTS.md). Its timestamped results live under
`reports/crc_experiments/`; the original compatibility outputs stay unchanged.

Run from the project root (paths are relative):

```sh
python3 -m pip install -r analysis/requirements.txt
python3 analysis/run_compatibility.py
python3 -m unittest discover -s tests/research -p 'test_compatibility.py'
```

The command inspects and hashes the eight requested ECoGraph exports and
`data/raw/primekg/kg.csv`, scans PrimeKG in chunks, audits identifiers and edges, computes
conservative symbol mappings, and generates `reports/compatibility/compatibility_report.md`.
It never merges graphs or modifies source datasets. It reuses the endpoint
schema and identifier-preservation approach of `src/primekg/prepare.py`, without
relying on possibly stale derived CSVs.

The supplied sample–gene export is blank. The default command additionally
verifies and analyzes the local `ECoGraph/data/coad_mut_data.csv` binary matrix
against exported sample, gene and subtype identities, checking subtype pairs
against `coad_subtype_data.csv`. It labels all resulting subtype gene sets and
metrics **supplementary**. It does not repair the export. Disable that additional
analysis with `--no-supplementary-recovery`; subtype molecular outputs will be
unavailable, not biologically empty. Use a different output directory for this
mode to avoid mixing outputs from different runs:

```sh
python3 analysis/run_compatibility.py --no-supplementary-recovery --output reports/compatibility_export_only
```

Other options: `--primekg`, `--ecograph`, and `--output`. No web services,
downloads, graph framework, statistical modeling package, or random seeds are
required. Results are deterministic for a given input. The primary raw graph is
not copied; only small audits and relevant direct target pairs are saved.

Mapping uses unique exact symbols to existing NCBI IDs, or whitespace trimming
with a separately recorded status. No synonym or case folding is accepted.
Pathway name candidates remain unverified; no fuzzy equivalences are accepted.

The supplementary statistical screens use one-sided Fisher/hypergeometric
tests and a one-sided, burden-quintile-stratified CMH normal approximation.
Each family receives BH correction over every subtype × eligible gene test.
Joint screen candidates are not validated enriched genes. The report explains
small-group, burden, and source limitations. Validated enrichment fields remain
unavailable. Missing metrics are written as empty CSV values, not zero.

Burden screens are disabled if either comparison group has fewer than 20
samples, or fewer than half of either group's samples lie in strata with both
groups represented. Disabled tests receive p=1 before BH; unfiltered approximate
p-values are diagnostic only. These thresholds are operational guardrails,
not proof that remaining asymptotic tests are reliable. In the supplied data,
both hypermutated groups fail these guardrails.

This workflow deliberately stops if the sample–gene export becomes populated:
its new schema and observation semantics must then be verified before replacing
the separately labeled recovery analysis. It does not silently reinterpret
new source data using assumptions from a blank export.

Start with `reports/compatibility/compatibility_report.md`; `reports/compatibility/README.md` lists artifacts.
Schemas/records are in `source_inspection.txt`; SHA-256/version evidence is in
`source_manifest.json`, `supplementary_recovery_audit.json`, and `validation.json`.
