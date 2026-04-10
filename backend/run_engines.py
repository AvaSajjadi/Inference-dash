from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def run_cie(signature_path: str, network_path: str, output_path: str, p_cut: str = "0.2", fc_cut: str = "1.1") -> None:
    signature_path = str(Path(signature_path))
    network_path = str(Path(network_path))
    output_path = str(Path(output_path))

    cmd = [
        "Rscript",
        "backend/engines/cie/run_cie.R",
        signature_path,
        network_path,
        output_path,
        str(p_cut),
        str(fc_cut),
    ]
    subprocess.run(cmd, check=True)


def run_ornor(signature_path: str, network_path: str, output_path: str, top_edges: int = 5000) -> None:
    signature_path = str(Path(signature_path))
    network_path = str(Path(network_path))
    output_path = str(Path(output_path))

    cmd = [
        sys.executable,
        "-m",
        "backend.engines.ornor.run_ornor_real",
        signature_path,
        network_path,
        output_path,
        str(int(top_edges)),
    ]
    subprocess.run(cmd, check=True)
