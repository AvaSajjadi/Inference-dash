from __future__ import annotations

import uuid
from pathlib import Path
from typing import Dict, Any

from backend.engines.base import BaseEngine
from backend.app.schema import (
    SignatureInput,
    NetworkInput,
    RunConfig,
    CanonicalOutputs,
)
from backend.app.provenance import write_provenance


class InferenceService:
    """
    High-level orchestration layer.
    Responsible for:
    - Creating run folder
    - Calling engine
    - Writing provenance
    """

    def __init__(self, results_root: Path):
        self.results_root = results_root

    def run(
        self,
        engine: BaseEngine,
        signature: SignatureInput,
        network: NetworkInput,
        config: RunConfig,
    ) -> CanonicalOutputs:

        run_id = str(uuid.uuid4())
        run_dir = self.results_root / "runs" / run_id
        run_dir.mkdir(parents=True, exist_ok=True)

        engine_output_dir = run_dir / "engine"
        engine_output_dir.mkdir(exist_ok=True)

        result = engine.run(
            signature_path=signature.path,
            network_path=network.path,
            output_dir=engine_output_dir,
            seed=config.seed,
            params=config.params,
        )

        provenance_path = run_dir / "provenance.json"

        write_provenance(
            output_path=provenance_path,
            signature_path=signature.path,
            network_path=network.path,
            engine_name=engine.name,
            engine_version=engine.version(),
            seed=config.seed,
            extra=result.diagnostics,
        )

        return CanonicalOutputs(
            edges_csv=result.edges_path,
            regulators_csv=result.regulators_path,
            diagnostics=result.diagnostics,
        )

