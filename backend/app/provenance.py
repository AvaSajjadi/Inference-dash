from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from pathlib import Path
from typing import Dict, Any
from datetime import datetime


def sha256_of_file(path: Path) -> str:
    """
    Compute SHA256 hash of a file for reproducibility tracking.
    """
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()


def get_python_version() -> str:
    return platform.python_version()


def get_r_version() -> str:
    try:
        result = subprocess.run(
            ["R", "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.split("\n")[0]
    except Exception:
        return "R not found"


def write_provenance(
    output_path: Path,
    signature_path: Path,
    network_path: Path,
    engine_name: str,
    engine_version: str,
    seed: int,
    extra: Dict[str, Any] | None = None,
):
    """
    Write full reproducibility metadata to JSON.
    """

    metadata = {
        "timestamp": datetime.utcnow().isoformat(),
        "engine": engine_name,
        "engine_version": engine_version,
        "seed": seed,
        "signature_sha256": sha256_of_file(signature_path),
        "network_sha256": sha256_of_file(network_path),
        "python_version": get_python_version(),
        "r_version": get_r_version(),
        "extra": extra or {},
    }

    with output_path.open("w") as f:
        json.dump(metadata, f, indent=2)
