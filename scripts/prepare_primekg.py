"""Run from the repository root: python3 scripts/prepare_primekg.py."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from primekg.prepare import main

if __name__ == "__main__":
    main()
