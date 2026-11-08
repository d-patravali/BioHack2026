"""Read only the Nebius settings; never expose the key in repr or logs."""

import os
from dataclasses import dataclass, field
from pathlib import Path

from dotenv import dotenv_values

PROJECT_ROOT = Path(__file__).resolve().parents[1]
NEBIUS_BASE_URL = "https://api.tokenfactory.nebius.com/v1/"


@dataclass(frozen=True)
class Settings:
    api_key: str = field(repr=False)
    model: str


def load_settings(env_file: Path = PROJECT_ROOT / ".env") -> Settings:
    values = {**dotenv_values(env_file), **os.environ}
    missing = [name for name in ("NEBIUS_API_KEY", "NEBIUS_MODEL")
               if not (values.get(name) or "").strip()]
    if missing:
        raise ValueError("Set " + ", ".join(missing) + " in the project .env or environment.")
    return Settings(api_key=values["NEBIUS_API_KEY"].strip(), model=values["NEBIUS_MODEL"].strip())
