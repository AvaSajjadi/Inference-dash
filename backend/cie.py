import subprocess
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

def run_cie(expr_csv_path):
    expr_csv_path = Path(expr_csv_path)

    network_csv = BASE_DIR / "networks" / "human_network.csv"
    out_csv = BASE_DIR / "results" / "cie_output.csv"

    cmd = [
        "Rscript",
        str(BASE_DIR / "run_cie.R"),
        str(expr_csv_path),
        str(network_csv),
        str(out_csv),
    ]

    subprocess.check_call(cmd)
    return out_csv
