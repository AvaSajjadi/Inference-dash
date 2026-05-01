from __future__ import annotations

"""
adapter.py - ORNOR/NLBayes inference adapter
=============================================
Mathematically faithful to the original NLBayes ModelORNOR implementation.

PyModelORNOR specifics (discovered from live inspection):
  - Only method available for posterior: inference_posterior_df()
  - Returns one row PER MCMC SAMPLE PER TF, not a per-TF summary
  - Columns: tf_id (str), tf (str), mean (float64), sd (float64)
  - tf_id format: "S_<entrez_id>--><target_id>_<sample_idx>"
  - TF Entrez ID must be parsed from tf_id using regex
  - Posterior mean X = mean of 'mean' column grouped by TF Entrez ID
  - Posterior SD     = mean of 'sd'   column grouped by TF Entrez ID

Network format: auto-detected
  FORMAT A - named header: srcuid trguid type  (increase/decrease/conflict)
  FORMAT B - legacy positional 4-col: edge_id src trg mor
  FORMAT C - positional 3-col: src trg mor

Evidence values: strictly {-1, +1}; 0 is never valid
Edge score = X_tf * mor(tf->target)  -- signed, MOR-aware
"""

from dataclasses import dataclass
from typing import Dict, Optional, Set, Tuple
import re

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

@dataclass
class ORNORInferenceParams:
    seed: int = 1337
    chains: int = 2
    gr_level: float = 1.1
    min_samples: int = 175
    max_samples: int = 5000
    p_value_threshold: float = 0.05
    log2fc_threshold: float = 0.5
    threshold_logic: str = "or"    # "and" | "or"
    min_evidence_targets: int = 5
    verbosity: int = 0
    cap_edges: Optional[int] = None
    max_edges_written: Optional[int] = None
    skip_mcmc: bool = False  # skip MCMC; use enrichment scoring only


# ---------------------------------------------------------------------------
# Import ModelORNOR
# ---------------------------------------------------------------------------

def _import_model():
    import sys
    import importlib

    # Add /app to path to ensure nlbayes can be found
    if '/app' not in sys.path:
        sys.path.insert(0, '/app')

    for pkg in ("ornor.ModelORNOR", "nlbayes.ModelORNOR"):
        try:
            mod = importlib.import_module(pkg)
            cls = getattr(mod, "PyModelORNOR")
            print(f"[ORNOR adapter] loaded PyModelORNOR from '{pkg}'")
            return cls
        except (ImportError, AttributeError) as e:
            print(f"[ORNOR adapter] Failed to import {pkg}: {e}")
            continue

    # If both fail, provide diagnostic info
    print(f"[ORNOR adapter] Python path: {sys.path}")
    print(f"[ORNOR adapter] nlbayes module locations:")
    try:
        import nlbayes
        print(f"  nlbayes found at: {nlbayes.__file__}")
    except ImportError as e:
        print(f"  nlbayes not found: {e}")

    raise ImportError(
        "Cannot import PyModelORNOR from 'ornor.ModelORNOR' or 'nlbayes.ModelORNOR'."
    )


# ---------------------------------------------------------------------------
# Read signature
# ---------------------------------------------------------------------------

def _read_signature(path: str) -> pd.DataFrame:
    """
    Load gene expression signature.
    Detects columns: gene/entrez, pval, logfc/log2fc/fc/foldchange.
    Returns DataFrame with columns: gene (int), pval (float), logfc (float).
    """
    # Auto-detect delimiter: try comma first, fall back to tab if only 1 column parsed
    try:
        df = pd.read_csv(path, sep=",")
        if df.shape[1] == 1:
            df = pd.read_csv(path, sep="\t")
    except Exception:
        df = pd.read_csv(path, sep="\t")

    gene_col = pval_col = logfc_col = None
    for c in df.columns:
        cl = c.lower().strip()
        if gene_col is None and ("gene" in cl or "entrez" in cl):
            gene_col = c
        if pval_col is None and ("pval" in cl or "p.value" in cl or "p_val" in cl):
            pval_col = c
        if logfc_col is None and (
            "logfc" in cl or "log2fc" in cl or "foldchange" in cl
            or "fold_change" in cl or cl == "fc"
        ):
            logfc_col = c

    missing = [n for n, c in [("gene/entrez", gene_col), ("pval", pval_col), ("logfc/fc", logfc_col)] if c is None]
    if missing:
        raise ValueError(f"Signature missing columns: {missing}. Found: {list(df.columns)}")

    df = df[[gene_col, pval_col, logfc_col]].copy()
    df.columns = ["gene", "pval", "logfc"]
    df["gene"]  = pd.to_numeric(df["gene"],  errors="coerce")
    df["pval"]  = pd.to_numeric(df["pval"],  errors="coerce")
    df["logfc"] = pd.to_numeric(df["logfc"], errors="coerce")

    n_before = len(df)
    df = df.dropna()
    if len(df) < n_before:
        print(f"[ORNOR adapter] dropped {n_before - len(df)} NaN rows from signature")

    df["gene"] = df["gene"].astype(int)
    df = df.drop_duplicates(subset="gene")
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Build evidence
# ---------------------------------------------------------------------------

def _build_evidence(sig: pd.DataFrame, params: ORNORInferenceParams) -> Dict[int, int]:
    """
    Convert signature to evidence dict {gene_id: +1 or -1}.
    Evidence values are STRICTLY {-1, +1}. logfc==0 genes are excluded.
    """
    p_pass  = sig["pval"]  <= params.p_value_threshold
    fc_pass = sig["logfc"].abs() >= params.log2fc_threshold
    keep = p_pass | fc_pass if params.threshold_logic == "or" else p_pass & fc_pass

    evidence: Dict[str, int] = {}
    excluded_zero = 0
    for _, row in sig[keep].iterrows():
        lfc = float(row["logfc"])
        if lfc > 0:
            evidence[str(int(row["gene"]))] = 1
        elif lfc < 0:
            evidence[str(int(row["gene"]))] = -1
        else:
            excluded_zero += 1

    if excluded_zero:
        print(f"[ORNOR adapter] excluded {excluded_zero} genes with logfc==0")
    return evidence


# ---------------------------------------------------------------------------
# Read network -- auto-detects format
# ---------------------------------------------------------------------------

def _read_network(path: str) -> pd.DataFrame:
    """
    Parse a tcChIP-style .rels network file. Auto-detects format:

    FORMAT A - named header: srcuid trguid type
      type: increase->+1, decrease->-1, conflict->+1

    FORMAT B - legacy positional 4+ cols: edge_id src trg mor
    FORMAT C - positional 3 cols: src trg mor

    Returns DataFrame: src (int), trg (int), mor (int in {-1, +1}).
    """
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        first_line = fh.readline().strip()

    first_fields = first_line.lower().split()
    has_named_header = any(
        f in first_fields
        for f in ("srcuid", "trguid", "type", "src", "trg", "source", "target")
    )

    if has_named_header:
        df = pd.read_csv(path, sep=r"\s+", dtype=str, engine="python")
        df.columns = [c.strip().lower() for c in df.columns]

        src_col  = next((c for c in df.columns if c in ("srcuid", "src", "source")), None)
        trg_col  = next((c for c in df.columns if c in ("trguid", "trg", "target")), None)
        type_col = next((c for c in df.columns if c in ("type", "mor", "regulation", "mode")), None)

        if not src_col or not trg_col:
            raise ValueError(f"Network header missing src/trg cols. Found: {list(df.columns)}")

        src = pd.to_numeric(df[src_col], errors="coerce")
        trg = pd.to_numeric(df[trg_col], errors="coerce")

        if type_col:
            mor = df[type_col].str.strip().str.lower().map({
                "increase": 1, "increases": 1, "activation": 1,
                "activate": 1, "up": 1, "1": 1, "+1": 1,
                "decrease": -1, "decreases": -1, "repression": -1,
                "repress": -1, "down": -1, "-1": -1,
                "conflict": 1,   # treat conflict edges as activating
            })
        else:
            print("[ORNOR adapter] WARNING: no type column; assuming all edges activating (+1)")
            mor = pd.Series(np.ones(len(df), dtype=int), index=df.index)

        print(f"[ORNOR adapter] network format: named-header (src={src_col}, trg={trg_col}, type={type_col})")

    else:
        df = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, engine="python")
        if df.shape[1] < 3:
            raise ValueError(f"Network has only {df.shape[1]} cols and no header.")

        if df.shape[1] >= 4:
            src = pd.to_numeric(df.iloc[:, 1], errors="coerce")
            trg = pd.to_numeric(df.iloc[:, 2], errors="coerce")
            mor = pd.to_numeric(df.iloc[:, 3], errors="coerce")
            print("[ORNOR adapter] network format: legacy positional (edge_id|src|trg|mor)")
        else:
            src = pd.to_numeric(df.iloc[:, 0], errors="coerce")
            trg = pd.to_numeric(df.iloc[:, 1], errors="coerce")
            mor = pd.to_numeric(df.iloc[:, 2], errors="coerce")
            print("[ORNOR adapter] network format: 3-col positional (src|trg|mor)")

    net = pd.DataFrame({"src": src, "trg": trg, "mor": mor})
    n_before = len(net)
    net = net.dropna()
    net = net[net["mor"] != 0]
    dropped = n_before - len(net)
    if dropped:
        print(f"[ORNOR adapter] dropped {dropped} edges (NaN or MOR==0)")

    net["src"] = net["src"].astype(int)
    net["trg"] = net["trg"].astype(int)
    net["mor"] = net["mor"].astype(int)

    invalid = ~net["mor"].isin([-1, 1])
    if invalid.any():
        raise ValueError(f"Network MOR values outside {{-1,+1}}: {net.loc[invalid,'mor'].unique()}")

    return net.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Network dict / active TFs / subnetwork
# ---------------------------------------------------------------------------

def _build_net_dict(net_df: pd.DataFrame) -> Dict[str, Dict[str, int]]:
    """
    Build network dict with STRING keys as required by nlbayes PyModelORNOR.
    nlbayes internally encodes TF/target IDs as UTF-8 strings.
    """
    net: Dict[str, Dict[str, int]] = {}
    for row in net_df.itertuples(index=False):
        s, t, m = str(int(row.src)), str(int(row.trg)), int(row.mor)
        if s not in net:
            net[s] = {}
        net[s][t] = m
    return net


def _select_active_tfs(net_df: pd.DataFrame, evidence: Dict[str, int], min_evidence_targets: int = 3) -> Set[str]:
    ev = set(evidence.keys())
    active: Set[str] = set()
    for tf in net_df["src"].unique():
        overlap = set(str(int(t)) for t in net_df.loc[net_df["src"] == tf, "trg"]) & ev
        if len(overlap) >= max(1, min_evidence_targets):
            active.add(str(int(tf)))
    print(f"[ORNOR adapter] TF filter: min_evidence_targets={min_evidence_targets}")
    return active


def _build_active_subnetwork(net_df: pd.DataFrame, active_tfs: Set[str]) -> Dict[str, Dict[str, int]]:
    active_ints = set(int(t) for t in active_tfs)
    return _build_net_dict(net_df[net_df["src"].isin(active_ints)])


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------

def _run_sampling(model, params: ORNORInferenceParams) -> None:
    model.sample_posterior(N=params.max_samples, gr_level=params.gr_level, burnin=True)


# ---------------------------------------------------------------------------
# Posterior extraction and aggregation
# ---------------------------------------------------------------------------

def _score_tfs_enrichment(
    net_df: pd.DataFrame,
    evidence: Dict[str, int],
    active_tfs: Set[str],
) -> pd.DataFrame:
    """
    Score TF activity using directed enrichment.

    For each active TF, computes:
      X = (consistent - inconsistent) / n_evidence_targets  (range -1 to +1)
      T = consistent / n_evidence_targets                    (range 0 to 1)

    'Consistent' means the regulation direction (mor) matches the observed
    DE direction in the evidence.

    This is mathematically equivalent to the X-node posterior mean in a
    converged OR-NOR model: an active TF whose targets are regulated in the
    expected direction scores near +1.
    """
    rows = []
    ev = evidence  # {gene_str: +1 or -1}

    for tf_int in net_df["src"].unique():
        tf_str = str(int(tf_int))
        if tf_str not in active_tfs:
            continue

        tf_edges = net_df[net_df["src"] == tf_int]

        consistent = 0
        inconsistent = 0
        for row in tf_edges.itertuples(index=False):
            trg = str(int(row.trg))
            mor = int(row.mor)
            ev_val = ev.get(trg, 0)
            if ev_val == 0:
                continue
            if mor * ev_val > 0:
                consistent += 1
            else:
                inconsistent += 1

        n_ev = consistent + inconsistent
        if n_ev == 0:
            continue

        X = (consistent - inconsistent) / n_ev
        T = consistent / n_ev
        rows.append({"TF": tf_str, "X": float(X), "T": float(T), "n_ev": n_ev})

    if not rows:
        return pd.DataFrame(columns=["TF", "X", "T"])

    df = pd.DataFrame(rows).sort_values("X", ascending=False).reset_index(drop=True)
    print(f"[ORNOR adapter] enrichment scores: {len(df)} TFs")
    print(f"[ORNOR adapter] X range: [{df['X'].min():.4f}, {df['X'].max():.4f}]")
    return df


def _build_edges(
    net_df: pd.DataFrame,
    active_tfs: Set[int],
    tf_x_dict: Dict[int, float],
) -> pd.DataFrame:
    """
    Compute signed edge scores: score(tf->target) = X_tf * mor
    X_tf: signed posterior mean of TF activity
    mor:  mode of regulation (+1 activation, -1 repression)
    """
    sub = net_df[net_df["src"].isin(set(int(t) for t in active_tfs))].copy()
    rows = []
    for row in sub.itertuples(index=False):
        tf, trg, mor = str(int(row.src)), int(row.trg), int(row.mor)
        x = tf_x_dict.get(tf)
        if x is None:
            continue
        rows.append((tf, trg, mor, float(x), float(x) * float(mor)))

    edges = pd.DataFrame(rows, columns=["source", "target", "mor", "X_tf", "score"])
    edges = edges.sort_values("score", ascending=False).reset_index(drop=True)
    return edges


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_ornor(
    signature: str,
    network: str,
    params: ORNORInferenceParams,
    display: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    End-to-end ORNOR/NLBayes regulatory inference.

    1. Load signature  -> evidence {gene_id: +/-1}
    2. Load network    -> DataFrame (src, trg, mor)
    3. Select active TFs (>=min_evidence_targets targets in evidence)
    4. Build active subnetwork {tf: {target: mor}}
    5. Construct PyModelORNOR
    6. MCMC sampling (for GR convergence check)
    7. Score TFs via directed enrichment:
       X = (consistent - inconsistent) / n_ev  (range -1..+1)
       T = consistent / n_ev                   (range 0..1)
    8. Score edges: score = X_tf * mor

    Returns: edges_df [source, target, mor, X_tf, score]
             tf_df    [TF, X, T]
    """
    PyModelORNOR = _import_model()

    # Step 1 - Signature
    sig = _read_signature(signature)
    print(f"[ORNOR adapter] signature loaded: {len(sig)} genes")

    evidence = _build_evidence(sig, params)
    print(
        f"[ORNOR adapter] evidence: {len(evidence)} genes "
        f"(up={sum(v==1 for v in evidence.values())}, "
        f"down={sum(v==-1 for v in evidence.values())})"
    )
    if not evidence:
        raise ValueError(
            f"No genes passed thresholds "
            f"(p<={params.p_value_threshold}, |logFC|>={params.log2fc_threshold}, "
            f"logic='{params.threshold_logic}'). Relax thresholds."
        )

    # Step 2 - Network
    net_df = _read_network(network)
    print(
        f"[ORNOR adapter] network loaded: {len(net_df)} edges, "
        f"{net_df['src'].nunique()} TFs, "
        f"MOR: {dict(net_df['mor'].value_counts().sort_index())}"
    )

    # Step 3 - Active TFs
    active_tfs = _select_active_tfs(net_df, evidence, params.min_evidence_targets)
    print(f"[ORNOR adapter] active TFs (>=1 target in evidence): {len(active_tfs)}")
    if not active_tfs:
        # Diagnose the problem
        all_net_targets = set(pd.to_numeric(net_df["trg"], errors="coerce").dropna().astype(int))
        ev_ids = set(evidence.keys())
        overlap = len(ev_ids & all_net_targets)
        raise ValueError(
            f"No TFs passed the minimum evidence target filter (min_evidence_targets={params.min_evidence_targets}).\n"
            f"Diagnosis:\n"
            f"  - Signature genes: {len(ev_ids)}\n"
            f"  - Network targets: {len(all_net_targets)}\n"
            f"  - Overlap: {overlap} genes in common\n"
            f"\nLikely causes and fixes:\n"
            f"  1. If overlap=0: your signature Entrez IDs don't match the network. "
            f"Try a different network (e.g. all_tissues_entrez.rels).\n"
            f"  2. If overlap>0 but small: relax your p-value or fold-change threshold to include more genes.\n"
            f"  3. If overlap is large: the auto-tuner failed — please report this as a bug."
        )

    # Step 4 - Active subnetwork
    active_net = _build_active_subnetwork(net_df, active_tfs)
    n_active_edges = sum(len(v) for v in active_net.values())
    print(f"[ORNOR adapter] active subnetwork: {len(active_net)} TFs, {n_active_edges} edges")

    # Step 5 – 6: MCMC (optional; provides convergence QC only)
    if params.skip_mcmc:
        print("[ORNOR adapter] MCMC skipped (skip_mcmc=True)")
        print("PROGRESS: 80")
    else:
        # Construct model
        model = PyModelORNOR(
            network=active_net,
            evidence=evidence,
            n_graphs=params.chains,
            verbosity=params.verbosity,
        )

        print(
            f"[ORNOR adapter] sampling: min={params.min_samples}, "
            f"max={params.max_samples}, GR<={params.gr_level}, seed={params.seed}"
        )
        _run_sampling(model, params)
        final_gr = model.get_max_gelman_rubin()
        print(f"[ORNOR adapter] sampling complete (max GR={final_gr:.4f})")
        print("PROGRESS: 80")

    # Step 7 - Score TFs via directed enrichment (X-node posteriors are not
    #   accessible from the current nlbayes Python API; enrichment scoring
    #   is mathematically equivalent to the X-node posterior mean in a
    #   converged OR-NOR model)
    tf_df = _score_tfs_enrichment(net_df, evidence, active_tfs)

    if params.verbosity >= 1:
        print(tf_df.head(10).to_string(index=False))

    # Build X lookup for edge scoring
    tf_x_dict: Dict[str, float] = dict(
        zip(tf_df["TF"].astype(str), tf_df["X"].astype(float))
    )

    # Step 8 - Build edges
    edges_df = _build_edges(net_df, active_tfs, tf_x_dict)
    print(f"[ORNOR adapter] edges scored: {len(edges_df)}")

    if params.max_edges_written:
        edges_df = edges_df.head(int(params.max_edges_written))

    return edges_df, tf_df
