# backend/ornor.py
"""
ORNOR stub backend.

This is a placeholder implementation that mimics the output
format of ORNOR inference without running nlbayes yet.
"""

from pathlib import Path
import pandas as pd


def run_ornor(expr_csv: Path, out_csv: Path) -> None:
    """
    Stub ORNOR inference.

    Parameters
    ----------
    expr_csv : Path
        Input expression matrix CSV
    out_csv : Path
        Output CSV path
    """

    # Fake ORNOR-style posterior results
    df = pd.DataFrame({
        "TF": [151636, 3799, 23317, 9652, 3895,
               55183, 23064, 84162, 2633, 7175],
        "posterior_mean": [0.91, 0.88, 0.86, 0.84, 0.82,
                           0.80, 0.78, 0.76, 0.74, 0.72],
        "posterior_sd": [0.02] * 10
    })

    df.to_csv(out_csv, index=False)
