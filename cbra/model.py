"""
CBRA: Continuous Bayesian Regulatory Activity model
for Toxoplasma gondii (and any organism with a custom network).

Model
-----
For each target gene g with observed log2FC y_g:

    y_g  ~  Normal( sum_k( a_k * W_kg ),  sigma )

    a_k  ~  Normal( 0,  tau * lambda_k )    [horseshoe prior]
    lambda_k ~ HalfCauchy(1)                [local shrinkage per TF]
    tau      ~ HalfCauchy(1)                [global sparsity]
    sigma    ~ HalfNormal(1)                [observation noise]

W_kg = +1  if TF k activates gene g
       -1  if TF k represses gene g
        0  no known edge

Output: posterior distribution over each TF's activity score a_k,
        ranked by |mean(a_k)| with 94% credible intervals.

Reference for horseshoe prior:
    Carvalho, Polson & Scott (2010). The horseshoe estimator for sparse signals.
    Biometrika, 97(2), 465-480.
"""

import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
import arviz as az
import matplotlib.pyplot as plt
from pathlib import Path


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_signature(sig_path: str | Path) -> pd.DataFrame:
    sig = pd.read_csv(sig_path, sep="\t")
    sig.columns = [c.lower() for c in sig.columns]

    id_col = next((c for c in sig.columns if c in ("entrez", "entrez_id", "gene_id", "geneid")), None)
    fc_col = next((c for c in sig.columns if c in ("fc", "logfc", "log2fc", "fold_change")), None)

    if id_col is None or fc_col is None:
        raise ValueError(f"Signature must have an ID column and a fold-change column. Found: {list(sig.columns)}")

    out = pd.DataFrame({
        "gene_id": sig[id_col].astype(str).str.strip(),
        "fc":      pd.to_numeric(sig[fc_col], errors="coerce"),
    }).dropna(subset=["fc"])

    out = out[out["gene_id"].ne("") & out["gene_id"].ne("nan")]
    return out.reset_index(drop=True)


def load_network(rels_path: str | Path, ents_path: str | Path):
    rels = pd.read_csv(rels_path, sep="\t")
    ents = pd.read_csv(ents_path, sep="\t")
    rels.columns = [c.lower() for c in rels.columns]
    ents.columns = [c.lower() for c in ents.columns]
    return rels, ents


# ---------------------------------------------------------------------------
# Weight matrix
# ---------------------------------------------------------------------------

def build_weight_matrix(
    sig: pd.DataFrame,
    rels: pd.DataFrame,
    ents: pd.DataFrame,
) -> tuple[np.ndarray, list[str], list[str], list[str]]:
    """
    Returns
    -------
    W          : (K, G) array  — +1 activation, -1 repression, 0 no edge
    tf_names   : list of K TF display names
    tf_ids     : list of K TF gene IDs
    target_ids : list of G target gene IDs (intersection of sig and network)
    """
    # Identify TFs (source nodes) and targets (sink nodes)
    src_uids = set(rels["srcuid"].unique())
    tf_ents  = ents[ents["uid"].isin(src_uids)].copy().reset_index(drop=True)

    # Map entity uid -> gene id string
    uid_to_id   = dict(zip(ents["uid"], ents["id"].astype(str)))
    uid_to_name = dict(zip(ents["uid"], ents["name"].astype(str)))

    # Relationship type -> sign
    def rel_sign(t):
        t = str(t).lower()
        if t in ("increase", "increases", "activation", "activates", "up"):
            return 1
        if t in ("decrease", "decreases", "repression", "represses", "down"):
            return -1
        return 0

    rels = rels.copy()
    rels["sign"] = rels["type"].apply(rel_sign)
    rels = rels[rels["sign"] != 0]

    # Target gene ids that appear in the signature
    sig_id_set   = set(sig["gene_id"])
    target_uids  = set(rels["trguid"].unique()) - src_uids
    target_ids_in_sig = [
        uid_to_id[u] for u in target_uids
        if uid_to_id.get(u, "") in sig_id_set
    ]
    target_ids_in_sig = sorted(set(target_ids_in_sig))

    if len(target_ids_in_sig) == 0:
        raise ValueError(
            "No overlap between network targets and signature gene IDs. "
            "Check that both files use the same ID format (e.g. TGME49_XXXXXX)."
        )

    gene_to_col = {gid: i for i, gid in enumerate(target_ids_in_sig)}
    K = len(tf_ents)
    G = len(target_ids_in_sig)
    W = np.zeros((K, G), dtype=np.float32)

    for _, row in rels.iterrows():
        src_uid = row["srcuid"]
        trg_id  = uid_to_id.get(row["trguid"], "")
        if trg_id not in gene_to_col:
            continue
        tf_k = tf_ents.index[tf_ents["uid"] == src_uid]
        if len(tf_k) == 0:
            continue
        W[tf_k[0], gene_to_col[trg_id]] = row["sign"]

    tf_names = [uid_to_name.get(u, str(u)) for u in tf_ents["uid"]]
    tf_ids   = [uid_to_id.get(u,   str(u)) for u in tf_ents["uid"]]

    print(f"[CBRA] {K} TF(s), {G} target genes in signature")
    print(f"[CBRA] W matrix density: {(W != 0).mean():.1%}")

    return W, tf_names, tf_ids, target_ids_in_sig


# ---------------------------------------------------------------------------
# PyMC model
# ---------------------------------------------------------------------------

def build_model(y: np.ndarray, W: np.ndarray) -> pm.Model:
    """
    y : (G,) observed log2FC for target genes
    W : (K, G) regulatory weight matrix
    """
    K, G = W.shape
    W_pt = pt.as_tensor_variable(W.astype(np.float32))

    with pm.Model() as model:
        # --- Horseshoe prior on TF activities ---
        # Global shrinkage: how sparse are active TFs overall?
        tau = pm.HalfCauchy("tau", beta=1.0)

        # Local shrinkage: allows a few TFs to escape global shrinkage
        lam = pm.HalfCauchy("lambda", beta=1.0, shape=K)

        # TF activity scores (continuous, signed)
        a = pm.Normal("a", mu=0.0, sigma=tau * lam, shape=K)

        # Predicted log2FC for each target gene: dot product of activities × edges
        mu = pt.dot(a, W_pt)   # shape (G,)

        # Observation noise (shared across genes; extend to per-gene if needed)
        sigma = pm.HalfNormal("sigma", sigma=1.0)

        # Likelihood
        pm.Normal("y_obs", mu=mu, sigma=sigma, observed=y.astype(np.float32))

    return model


# ---------------------------------------------------------------------------
# Inference
# ---------------------------------------------------------------------------

def run_inference(
    model: pm.Model,
    draws: int = 2000,
    tune:  int = 1000,
    chains: int = 2,
    target_accept: float = 0.9,
) -> az.InferenceData:
    with model:
        trace = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            target_accept=target_accept,
            progressbar=True,
            return_inferencedata=True,
        )
    return trace


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

def summarize(
    trace: az.InferenceData,
    tf_names: list[str],
    tf_ids: list[str],
) -> pd.DataFrame:
    summary = az.summary(trace, var_names=["a"], hdi_prob=0.94)
    summary.index = tf_names

    summary.insert(0, "tf_id", tf_ids)
    summary.insert(0, "tf_name", tf_names)
    summary = summary.rename(columns={
        "mean": "activity_mean",
        "sd":   "activity_sd",
        "hdi_3%":  "hdi_low",
        "hdi_97%": "hdi_high",
    })
    summary["abs_activity"] = summary["activity_mean"].abs()
    summary = summary.sort_values("abs_activity", ascending=False)
    return summary


def plot_results(summary: pd.DataFrame, out_path: str | Path | None = None):
    fig, ax = plt.subplots(figsize=(7, max(3, len(summary) * 0.5)))

    names  = summary["tf_name"]
    means  = summary["activity_mean"]
    lo     = summary["hdi_low"]
    hi     = summary["hdi_high"]
    colors = ["#e05252" if m > 0 else "#5278e0" for m in means]

    y_pos = range(len(names))
    ax.barh(y_pos, means, color=colors, alpha=0.75, height=0.5)
    ax.errorbar(
        means, y_pos,
        xerr=[means - lo, hi - means],
        fmt="none", color="black", capsize=4, linewidth=1.2,
    )
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(names)
    ax.set_xlabel("Activity score (posterior mean ± 94% HDI)")
    ax.set_title("CBRA: TF Regulatory Activity")
    plt.tight_layout()

    if out_path:
        fig.savefig(out_path, dpi=150)
        print(f"[CBRA] Plot saved to {out_path}")
    else:
        plt.show()

    return fig


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_cbra(
    sig_path:  str | Path,
    rels_path: str | Path,
    ents_path: str | Path,
    draws:   int = 2000,
    tune:    int = 1000,
    chains:  int = 2,
    out_dir: str | Path | None = None,
):
    sig        = load_signature(sig_path)
    rels, ents = load_network(rels_path, ents_path)
    W, tf_names, tf_ids, target_ids = build_weight_matrix(sig, rels, ents)

    # Align y to the target genes in W
    sig_indexed = sig.set_index("gene_id")
    y = np.array([sig_indexed.loc[gid, "fc"] for gid in target_ids], dtype=np.float32)

    print(f"[CBRA] y range: [{y.min():.2f}, {y.max():.2f}]  mean={y.mean():.2f}")

    model = build_model(y, W)
    trace = run_inference(model, draws=draws, tune=tune, chains=chains)
    summary = summarize(trace, tf_names, tf_ids)

    print("\n--- CBRA Results ---")
    print(summary[["tf_name", "tf_id", "activity_mean", "activity_sd", "hdi_low", "hdi_high"]].to_string())

    if out_dir:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out_dir / "cbra_results.tsv", sep="\t", index=False)
        trace.to_netcdf(str(out_dir / "cbra_trace.nc"))
        plot_results(summary, out_dir / "cbra_activity.png")
        print(f"[CBRA] Results written to {out_dir}")

    return trace, summary


if __name__ == "__main__":
    SIG   = Path("~/Downloads/toxo_bradyzoite_vs_tachyzoite_sig.tsv").expanduser()
    RELS  = Path("~/Downloads/toxo networks/toxo_ap2xii9.rels").expanduser()
    ENTS  = Path("~/Downloads/toxo networks/toxo_ap2xii9.ents").expanduser()
    OUT   = Path("~/inference_dash/cbra/results").expanduser()

    trace, summary = run_cbra(SIG, RELS, ENTS, draws=2000, tune=1000, chains=2, out_dir=OUT)
