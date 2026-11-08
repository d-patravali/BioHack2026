# Offline automated tests

Run from the repository root:

```sh
# biohack-agent environment:
python -m unittest discover -s tests/agent -p 'test_*.py' -v
# Research environment with pandas, NumPy, SciPy:
python3 -m unittest discover -s tests/research -p 'test_*.py' -v
```

An environment with all root requirements and pytest can use `python -m pytest`.
`pytest.ini` restricts discovery to this directory. None of these tests invokes
a live API, GNN training, or graph merging. Explicit live API demonstrations are
in `scripts/smoke/`.

CRC diagnostic tests also use optional PyTorch/scikit-learn packages and skip
when those packages are absent. Their tiny synthetic fixtures do not train a
GNN or access large research datasets.
