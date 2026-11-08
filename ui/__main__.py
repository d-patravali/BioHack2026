"""Start the local demonstration: python -m ui [--port 8765]."""

import argparse
import uvicorn

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    uvicorn.run("ui.server:app", host="127.0.0.1", port=args.port, access_log=False)
