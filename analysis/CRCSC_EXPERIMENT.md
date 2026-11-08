# CRCSC expression-profile reproducibility experiment

This is a separate CMS1–CMS4 representation experiment on `CRCSC-branch`.
It preserves the existing ECoGraph/CRC agent and the breast-cancer notebook.
It does **not** fit a drug predictor or GNN, merge a full graph, or infer treatments.

## Inputs

User-supplied files (read only, ignored by Git):

- `data/CRCSC_data/GSE39582_frma_expression.tsv`
- `data/CRCSC_data/FRENCH_clinical.tsv`
- `data/CRCSC_data/cms_labels_public_all.txt`

Other read-only inputs:

- `data/processed/primekg/nodes.csv`
- `data/raw/primekg/kg.csv`
- `data/external/geo/GPL570.annot.gz`

The official annotation was downloaded from
https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz
and its gzip integrity verified. Metadata dates it to Aug 09 2016. It is a
historical annotation snapshot; ambiguous probes are excluded rather than
silently assigned to one gene. Source hashes are recorded per completed run.
The script never downloads files automatically.

## Run

Use the research environment with pandas, NumPy, SciPy, scikit-learn, and
threadpoolctl. PyTorch and API credentials are not required.

```sh
python3 -m unittest discover -s tests/research -p 'test_crcsc_experiment.py' -v
python3 analysis/crcsc_experiment.py
python3 analysis/crcsc_report.py reports/crcsc/RUN_ID
```

A new UTC timestamped directory is created under `reports/crcsc/`. A supplied
`--output` must not exist. `--annotation` and `--input` override input locations.
The prespecified signature size is 100 (`--signature-size`). All models and
feature selection are fitted on labeled discovery samples. Validation labels
are used only to report results. Repeated tuning against these results would
require a fresh test cohort.

## Methods and boundaries

The expression header has no probe-ID field name, while every data row includes
that extra field. The loader handles this explicitly. GSM identifiers join the
three files one-to-one. Final combined CMS labels are used; NOLBL is not assigned
to a class. Multiple unambiguous probes per NCBI gene are aggregated by fixed
within-sample median. Genes map to PrimeKG by exact NCBI ID, not fuzzy names.

Discovery-standardized expression yields signed subtype-vs-rest signatures,
with OLS adjustment for batch, age, gender, stage, and tumor location. Discovery
BH q<=0.05 and absolute adjusted effect determine up to 100 selected genes per
subtype. Validation effect correlation, direction, and adjusted statistical
support assess replication. These are not tumor-versus-normal contrasts or
therapeutic dependencies.

Pathway scores average expression over measured Reactome members (>=5 genes).
Logistic and nearest-centroid classifiers assess CMS reproduction. Compare with
majority/batch baselines, expression-only logistic regression, and three gene-
permuted pathway controls. These controls preserve pathway sizes and membership
structure but do not preserve the degree of each named gene. Bootstrap intervals
for accuracy differences resample validation patients within CMS classes; they
cannot resolve circularity of expression-derived CMS labels or dataset bias.

CMS label reproduction is not independent subtype discovery, and no GNN
superiority is evaluated here. Drug-target connectivity is descriptive only.
The existing broad-disease indication task and shared-target drug ties are not
fixed by changing the molecular input. Downstream subtype–drug inference still
needs suitable supervision, drug distinctions, and independent response labels.
Literature review can check support, contradiction, and prior knowledge, not
independently verify efficacy.

## Outputs

Read `assessment.md`, `findings.md`, and `manifest.json` in the completed run first. Small CSVs
capture counts, joins/mappings, selected signatures, pathway profiles,
validation predictions, classifier metrics, paired bootstrap differences, and
direct-target connectivity. Large profile/audit tables are gzip-compressed.
No full processed expression matrix is duplicated on disk. A run without a
completed manifest is not a finished result.
