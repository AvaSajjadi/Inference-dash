from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, Optional


# ================================
# Canonical Input Definitions
# ================================

@dataclass
class SignatureInput:
    """
    Canonical representation of a signature file.
    Must already be validated and normalized.
    """
    path: Path


@dataclass
class NetworkInput:
    """
    Canonical representation of a network file.
    """
    path: Path


@dataclass
class RunConfig:
    """
    Configuration for a single inference run.
    """
    method: str  # "CIE", "ORNOR", or "BOTH"
    seed: int
    params: Dict[str, Any]


# ================================
# Canonical Output Definitions
# ================================

@dataclass
class CanonicalOutputs:
    """
    Standardized output files that UI will consume.
    """
    edges_csv: Path
    regulators_csv: Optional[Path]
    diagnostics: Dict[str, Any]
