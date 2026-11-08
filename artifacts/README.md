# Model artifacts

Generated files belong here rather than beside source code or raw datasets:

- `predictions/colorectal/`: CRC ranking exports from the notebook.
- `predictions/breast/`: breast-cancer mapping audit, predictions and rankings.
- `checkpoints/`: reserved for future saved model weights; no saving logic added.

These generated directories are ignored by Git. No models were trained during
the reorganization. Graph CSVs instead belong under `data/processed/`.
