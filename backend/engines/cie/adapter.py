# backend/engines/ornor/adapter.py
from __future__ import annotations

import os
import random
from dataclasses import dataclass
from typing import Optional, Tuple, Dict

import numpy as np
import pandas as pd


# -----------------------------
# Params
# -----------------------------
@dataclass
class ORNORInferenceParams:
    seed: int = 1337
    chains: int = 1
    gr_level: float = 1.1
    min_samples: int = 500
    max_samples: int = 5000

    # thresholds
    p_value_threshold: float = 0.05
    log2fc_threshold: float = 0.5
    threshold_logic: str = "and"  # used only if p-values exist

    professor_style: bool = True
    enable_fallback: bool = True

    cap_edges: Optional[int] = None
    max_edges_written: int = 5000


@dataclass
class ORNORDisplayParams:
    top_n_tfs_preview: int = 50
    top_n_edges_preview: int = 200
    min_tf_posterior_preview: float = 0.0


# -----------------------------
# Utilities
# -----------------------------
def _read_table(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path)
        if df.shape[1] == 1:
            df = pd.read_csv(path, sep="\t")
        return df
    except Exception:
        return pd.read_csv(path, sep="\t")


def _set_seeds(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def _infer_signature_columns(sig: pd.DataFrame) -> Dict[str, Optional[str]]:
    sig = sig.rename(columns={c: c.strip().lower() for c in sig.columns})
    cols = list(sig.columns)

    def pick(names):
        for n in names:
            if n in cols:
                return n
        return None

    gene_col = pick(["entrez", "gene", "geneid", "id", "uid", "ncbi", "ncbi_id"])
    pval_col = pick(["pval", "pvalue", "p.value", "padj", "adj.pval", "fdr", "qvalue", "e-fdr"])
    effect_col = pick(["fc", "log2fc", "log2foldchange", "log2fdc", "log2_fdc", "effect"])
    state_col = pick(["direction", "state", "sign"])

    if gene_col is None and cols:
        gene_col = cols[0]

    return {"gene_col": gene_col, "pval_col": pval_col, "effect_col": effect_col, "state_col": state_col}


def _read_network_rels(network_path: str) -> pd.DataFrame:
    rels = _read_table(network_path)
    rels = rels.rename(columns={c: c.strip().lower() for c in rels.columns})

    source_candidates = [c for c in rels.columns if c in ("source", "src", "tf", "regulator")]
    target_candidates = [c for c in rels.columns if c in ("target", "tgt", "gene", "targetgene")]

    if not source_candidates or not target_candidates:
        if rels.shape[1] >= 2:
            rels = rels.iloc[:, :3].copy()
            rels.columns = ["source", "target"] + (["mode"] if rels.shape[1] >= 3 else [])
        else:
            raise ValueError(f"Network rels missing source/target columns: {list(rels.columns)}")

    if "source" not in rels.columns:
        rels = rels.rename(columns={source_candidates[0]: "source"})
    if "target" not in rels.columns:
        rels = rels.rename(columns={target_candidates[0]: "target"})

    keep = ["source", "target"] + (["mode"] if "mode" in rels.columns else [])
    rels = rels[keep].copy()
    rels["source"] = rels["source"].astype(str)
    rels["target"] = rels["target"].astype(str)
    return rels


def _build_evidence(
    sig: pd.DataFrame,
    gene_col: str,
    pval_col: Optional[str],
    effect_col: str,
    params: ORNORInferenceParams,
) -> pd.DataFrame:
    """
    If pval_col is None -> filter only by FC threshold.
    If pval_col exists -> combine with threshold_logic.
    """
    fc = pd.to_numeric(sig[effect_col], errors="coerce").fillna(0.0)
    keep_fc = fc.abs() >= float(params.log2fc_threshold)

    if pval_col is None:
        keep = keep_fc
    else:
        pv = pd.to_numeric(sig[pval_col], errors="coerce").fillna(1.0)
        keep_pv = pv <= float(params.p_value_threshold)
        logic = (params.threshold_logic or "and").lower().strip()
        keep = (keep_fc | keep_pv) if logic == "or" else (keep_fc & keep_pv)

    direction = fc.apply(lambda x: 1 if x > 0 else (-1 if x < 0 else 0))
    direction = direction.where(keep, 0).astype(int)

    ev = pd.DataFrame({"gene": sig[gene_col].astype(str), "direction": direction})
    return ev.loc[ev["direction"] != 0].copy()


def _overlap_fallback(
    rels: pd.DataFrame,
    evidence_map: Dict[str, int],
    out_edges_csv: str,
    out_tf_tsv: str,
    max_edges: int,
) -> Tuple[str, str]:
    """
    Fallback TF ranking: overlap count of evidence genes with each TF's targets.
    Writes out_tf_tsv and out_edges_csv, then returns paths.
    """
    ev_genes = set(map(str, evidence_map.keys()))
    grp = rels.groupby("source")["target"].apply(list)

    tf_rows = []
    for tf, targets in grp.items():
        overlap = len(ev_genes.intersection(set(map(str, targets))))
        tf_rows.append((str(tf), float(overlap)))

    tf_post = (
        pd.DataFrame(tf_rows, columns=["TF", "posterior_mean"])
        .sort_values(["posterior_mean", "TF"], ascending=[False, True])
        .reset_index(drop=True)
    )

    os.makedirs(os.path.dirname(out_tf_tsv) or ".", exist_ok=True)
    tf_post.to_csv(out_tf_tsv, sep="\t", index=False)

    edges = []
    for tf, score in tf_post[["TF", "posterior_mean"]].itertuples(index=False, name=None):
        sub = rels.loc[rels["source"].astype(str) == str(tf)]
        for _, r in sub.iterrows():
            mode_val = int(r["mode"]) if "mode" in sub.columns else 0
            edges.append((str(r["source"]), str(r["target"]), float(score), mode_val))
            if len(edges) >= max_edges:
                break
        if len(edges) >= max_edges:
            break

    edges_df = pd.DataFrame(edges, columns=["source", "target", "score", "mode"])
    os.makedirs(os.path.dirname(out_edges_csv) or ".", exist_ok=True)
    edges_df.to_csv(out_edges_csv, index=False)

    print(f"[ORNOR adapter] FALLBACK wrote TFs : {out_tf_tsv} (rows={len(tf_post)})")
    print(f"[ORNOR adapter] FALLBACK wrote edges: {out_edges_csv} (rows={len(edges_df)})")
    return out_edges_csv, out_tf_tsv


# -----------------------------
# Main entry
# -----------------------------
def run_ornor(
    signature_path: str,
    network_path: str,
    out_edges_csv: str,
    out_tf_tsv: str,
    params: ORNORInferenceParams,
    display: Optional[ORNORDisplayParams] = None,
) -> Tuple[str, str]:
    _set_seeds(int(params.seed))
    if display is None:
        display = ORNORDisplayParams()

    # Read inputs
    sig_raw = _read_table(signature_path)
    rels = _read_network_rels(network_path)

    cols = _infer_signature_columns(sig_raw)
    gene_col = cols["gene_col"]
    pval_col = cols["pval_col"]
    effect_col = cols["effect_col"]

    if gene_col is None or effect_col is None:
        raise ValueError(f"Signature must have gene + fc/log2fc columns. Found: {list(sig_raw.columns)}")

    sig = sig_raw.copy()
    sig.columns = [c.strip().lower() for c in sig.columns]

    gene_col = gene_col.lower()
    if pval_col is not None:
        pval_col = pval_col.lower()
    effect_col = effect_col.lower()

    # Build evidence
    evidence = _build_evidence(sig, gene_col, pval_col, effect_col, params)
    evidence_map: Dict[str, int] = dict(zip(evidence["gene"].astype(str), evidence["direction"].astype(int)))

    # Logging (matches what you see)
    print(f"[ORNOR adapter] signature: {signature_path}")
    print(f"[ORNOR adapter] network: {network_path}")
    print(f"[ORNOR adapter] out_edges: {out_edges_csv}")
    print(f"[ORNOR adapter] out_tf: {out_tf_tsv}")
    print(f"[ORNOR adapter] inference params: {params}")
    print(f"[ORNOR adapter] display params (UI-only): {display}")
    print(f"[ORNOR adapter] signature columns detected: {cols}")
    print(f"[ORNOR adapter] signature genes (rows): {len(sig)}")
    print(f"[ORNOR adapter] evidence entries (nonzero): {len(evidence)}")
    print(f"[ORNOR adapter] TFs in network: {rels['source'].nunique()}")
    print(f"[ORNOR adapter] edges in network: {len(rels)}")

    # If no evidence -> fallback (or error if disabled)
    if len(evidence_map) == 0:
        if not params.enable_fallback:
            raise RuntimeError("No evidence after thresholding and fallback disabled.")
        print("[ORNOR adapter] No evidence -> fallback overlap ranking.")
        return _overlap_fallback(
            rels=rels,
            evidence_map={},
            out_edges_csv=out_edges_csv,
            out_tf_tsv=out_tf_tsv,
            max_edges=int(params.max_edges_written),
        )

    # Import nlbayes
    try:
        from nlbayes.ModelORNOR import PyModelORNOR  # type: ignore
    except Exception as e:
        print(f"[ORNOR adapter] Could not import PyModelORNOR: {type(e).__name__}: {e}")
        if not params.enable_fallback:
            raise
        print("[ORNOR adapter] Import failed -> fallback overlap ranking.")
        return _overlap_fallback(
            rels=rels,
            evidence_map=evidence_map,
            out_edges_csv=out_edges_csv,
            out_tf_tsv=out_tf_tsv,
            max_edges=int(params.max_edges_written),
        )

    # Build network dataframe expected by wrapper
    network_df = rels.copy()
    if "mode" not in network_df.columns:
        network_df["mode"] = 0

    # ✅ CRITICAL FIX: catch constructor/sampling/posterior crash and fallback instead of dying
    try:
        model = PyModelORNOR(
            network=network_df,
            evidence=evidence_map,
        )

        # Sampling
        model.set_seed(int(params.seed))
        model.sample_until_converged(
            gr_level=float(params.gr_level),
            min_samples=int(params.min_samples),
            max_samples=int(params.max_samples),
            chains=int(params.chains),
        )

        # Posterior TF table
        tf_post = model.inference_posterior_df()  # expected DataFrame

        # Write TF table to TSV
        os.makedirs(os.path.dirname(out_tf_tsv) or ".", exist_ok=True)
        tf_post.to_csv(out_tf_tsv, sep="\t", index=False)

        # Edge scoring
        try:
            edges_df = model.edge_scores_df()  # type: ignore
        except Exception:
            tf_col = "TF" if "TF" in tf_post.columns else ("tf" if "tf" in tf_post.columns else None)
            score_col = None
            for cand in ["posterior_mean", "mean", "prob", "score"]:
                if cand in tf_post.columns:
                    score_col = cand
                    break

            if tf_col is None or score_col is None:
                raise RuntimeError("Missing TF/score columns in posterior table for edge scoring.")

            tf_score = dict(
                zip(
                    tf_post[tf_col].astype(str),
                    pd.to_numeric(tf_post[score_col], errors="coerce").fillna(0.0),
                )
            )
            edges_df = rels.copy()
            edges_df["score"] = edges_df["source"].map(tf_score).fillna(0.0)
            if "mode" not in edges_df.columns:
                edges_df["mode"] = 0
            edges_df = edges_df[["source", "target", "score", "mode"]]

        if len(edges_df) > int(params.max_edges_written):
            edges_df = edges_df.head(int(params.max_edges_written)).copy()

        os.makedirs(os.path.dirname(out_edges_csv) or ".", exist_ok=True)
        edges_df.to_csv(out_edges_csv, index=False)

        print(f"[ORNOR adapter] wrote TFs : {out_tf_tsv} (rows={len(tf_post)})")
        print(f"[ORNOR adapter] wrote edges: {out_edges_csv} (rows={len(edges_df)})")
        return out_edges_csv, out_tf_tsv

    except Exception as e:
        print(f"[ORNOR adapter] NLBayes crashed during init/sampling/posterior: {type(e).__name__}: {e}")
        if not params.enable_fallback:
            raise
        print("[ORNOR adapter] Falling back to overlap ranking (stable output).")
        return _overlap_fallback(
            rels=rels,
            evidence_map=evidence_map,
            out_edges_csv=out_edges_csv,
            out_tf_tsv=out_tf_tsv,
            max_edges=int(params.max_edges_written),
        )
