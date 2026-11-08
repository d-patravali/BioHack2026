# Data layout

- `raw/primekg/kg.csv`: unchanged original PrimeKG input; not committed.
- `processed/primekg/`: initial node and phenotype–disease–gene exports.
- `processed/primekg_ecograph/`: reserved for the integration notebook's graph.
- `external/hgnc/`: versioned reference downloads used by integration.
- `external/gdsc/`: GDSC workbook used by the breast-cancer experiment.

Large raw, external, and processed data are ignored by Git. No graph merge was
executed as part of repository organization. The upstream ECoGraph data remain
in the separate root `ECoGraph/` checkout. See
[migration notes](../docs/repository_structure.md) for historical paths.
