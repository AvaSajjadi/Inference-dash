"""
ORNOR backend wrapper for Dash
- Ensures nlbayes-python is importable
- Runs ORNOR inference
- Returns posterior edge probabilities in a CIE-compatible format
"""

import sys
from pathlib import Path
import pandas as pd

# ------------------------------------------------------------------
# FIX PATH so Dash can import nlbayes even outside nlbayes-python
# ------------------------------------------------------------------

NLBAYES_ROOT = Path("/home/ava/nlbayes_project/nlbayes-python")
if str(NLBAYES_ROOT) not in sys.path:
    sys.path.insert(0, str(NLBAYES_ROOT))

from nlbayes import ORNOR


# ------------------------------------------------------------------
# ORNOR runner
# ------------------------------------------------------------------

def run_ornor(
    expression_csv: str,
    out_dir: str,
    n_samples: int = 5,
    burnin: int = 20,
):
    """
    Run ORNOR Bayesian inference

    Parameters
    ----------
    expression_csv : str
        Path to expression matrix CSV
    out_dir : str
        Output directory
    n_samples : int
        Posterior samples
    burnin : int
        Burn-in samples

    Returns
    -------
    pandas.DataFrame
        Posterior edge probabilities (CIE-style)
    """

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ----------------------------
    # Load expression data
    # ----------------------------
    expr = pd.read_csv(expression_csv, index_col=0)

    # ----------------------------
    # Initialize ORNOR
    # ----------------------------
    model = ORNOR(expr)

    # ----------------------------
    # Run inference
    # ----------------------------
    model.fit(expr, n_samples=n_samples, burnin=burnin)

    # ----------------------------
    # Extract posterior edges
    # ----------------------------
    edges = model.get_posterior()

    # Convert to DataFrame
    df = pd.DataFrame(
        edges,
        columns=["regulator", "target", "state", "probability"]
    )

    # Save CSV
    out_csv = out_dir / "ornor_posterior.csv"
    df.to_csv(out_csv, index=False)

    return df
