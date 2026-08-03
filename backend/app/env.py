"""Load `.env` at startup.

Kept dependency free on purpose: the server env is deliberately thin, and a
five line parser is less to justify than another package. `.env` is gitignored;
`.env.example` shows the shape.

Real environment variables always win, so an exported key overrides the file.
"""

from __future__ import annotations

import os
from pathlib import Path

ENV_FILE = Path(__file__).resolve().parents[2] / ".env"


def load(path: Path = ENV_FILE) -> list[str]:
    """Read KEY=value lines into os.environ. Returns the names that were set."""
    if not path.is_file():
        return []
    loaded: list[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        key, value = key.strip(), value.strip().strip("\'\"")
        # An exported variable beats the file, so a shell override still works.
        if key and key not in os.environ:
            os.environ[key] = value
            loaded.append(key)
    return loaded
