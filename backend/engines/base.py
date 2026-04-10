from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Any


class EngineResult:
    """
    Canonical engine result container.
    This is what ALL engines must return.
    """

    def __init__(
        self,
        edges_path: Path,
        regulators_path: Path | None,
        diagnostics: Dict[str, Any],
        raw_output_dir: Path,
    ):
        self.edges_path = edges_path
        self.regulators_path = regulators_path
        self.diagnostics = diagnostics
        self.raw_output_dir = raw_output_dir


class BaseEngine(ABC):
    """
    Abstract base class for all inference engines.
    Ensures deterministic and reproducible execution.
    """

    name: str

    @abstractmethod
    def run(
        self,
        signature_path: Path,
        network_path: Path,
        output_dir: Path,
        seed: int,
        params: Dict[str, Any],
    ) -> EngineResult:
        """
        Execute the engine.

        Must:
        - NOT reinterpret math
        - Call the real engine backend
        - Write raw outputs
        - Produce canonical edges CSV
        - Return EngineResult
        """
        raise NotImplementedError

    @abstractmethod
    def version(self) -> str:
        """
        Return version string for reproducibility logging.
        """
        raise NotImplementedError
