# Gene–phenotype evidence robustness experiment

This is separate from the working CRC agent/demo and preserves existing inputs and results.

```sh
python3 -m unittest discover -s tests/research -p 'test_phenotype_gnn.py' -v
python3 analysis/phenotype_gnn.py
python3 analysis/phenotype_report.py reports/phenotype_gnn/RUN_ID
```

Use a research environment with pandas, NumPy, SciPy, and PyTorch. The minimal `biohack-agent` environment is not intended for training. CPU execution uses two threads; no external service or graph library is needed. Defaults: 60 epochs, width24, masks17/29/43, initialization42. `--output` must not exist; `--epochs`, `--width`, `--seeds`, and `--evaluate-every` support debugging. Source files stay read only. A new timestamped results folder is created under `reports/phenotype_gnn/`.

Read that run's `findings.md` for actual results. It compares a two-layer relational GNN against matched popularity, PPI, pathway, fixed fusion, and disease-path heuristics. Entire annotation sets are withheld for query genes; validation genes are separate. Existing historically inspected test masks are reused, so results are exploratory, not a fresh confirmatory benchmark. No disease-path-derived labels are created.

The experiment also removes evidence channels at inference and evaluates validation-selected confidence-margin thresholds. These diagnostics establish neither causal explanations nor biological truth. The runner does not modify the agent, register a new model tool, merge the full graph, or claim a GNN advantage. Graph labels are sparse and missing edges remain unknown.

The report generator requires a completed preservation audit. It can regenerate narrative/summary tables without retraining. Saved checkpoints contain only tensors and simple metadata; load only trusted local artifacts.
