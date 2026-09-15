"""Path safety for the data-load endpoints.

The Data tab sends a filesystem path typed by the user. We only ever read files
under ``$SJ_DATA_ROOT`` (default: the user's home directory) and never follow a
path that escapes it.
"""
from __future__ import annotations

import os
from pathlib import Path

from fastapi import HTTPException


def data_root() -> Path:
    return Path(os.environ.get("SJ_DATA_ROOT", str(Path.home()))).expanduser().resolve()


def safe_path(raw: str) -> Path:
    if not raw or not raw.strip():
        raise HTTPException(400, "a path is required")
    try:
        p = Path(raw.strip()).expanduser().resolve()
    except OSError:
        raise HTTPException(400, f"cannot resolve path {raw!r}")
    root = data_root()
    if p != root and root not in p.parents:
        raise HTTPException(403, f"path must be inside {root}")
    if not p.exists():
        raise HTTPException(404, f"no such path: {p}")
    return p
