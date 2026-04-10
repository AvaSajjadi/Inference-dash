from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Dict

import pandas as pd

from backend.engines.ornor.adapter import ORNORInferenceParams, run_ornor


def _load_entrez_symbol_map(entities_path: str) -> Dict[str, str]:
    """Load uid->name map from a tcChIP .entities file (tab-separated, uid/name cols)."""
    try:
        df = pd.read_csv(entities_path, sep="\t", dtype=str)
        df.columns = [c.strip().lower() for c in df.columns]
        uid_col = next((c for c in df.columns if c in ("uid", "id", "entrez", "entrezid")), None)
        name_col = next((c for c in df.columns if c in ("name", "symbol", "gene_symbol", "label")), None)
        if not uid_col or not name_col:
            print(f"[ORNOR runner] entities: could not find uid/name cols in {entities_path} (found: {list(df.columns)})")
            return {}
        mapping = dict(zip(df[uid_col].str.strip(), df[name_col].str.strip()))
        print(f"[ORNOR runner] loaded {len(mapping)} symbol mappings from {entities_path}")
        return mapping
    except Exception as e:
        print(f"[ORNOR runner] WARNING: failed to load entities file {entities_path}: {e}")
        return {}


def main() -> int:
    ap = argparse.ArgumentParser(description="Run ORNOR (NLBayes) inference end-to-end.")
    ap.add_argument("signature", help="Signature file (tsv/csv): must include entrez,pval,fc")
    ap.add_argument("network", help="Network file (.rels): expected tcChIP style (edge_id,src,trg,mor,...)")
    ap.add_argument("out_edges", help="Output edges CSV path")
    ap.add_argument("out_tfs", help="Output TF TSV path")
    ap.add_argument("top_edges", type=int, help="Number of top edges to keep")

    ap.add_argument("--seed", type=int, default=1337)
    ap.add_argument("--chains", type=int, default=2)
    ap.add_argument("--gr_level", type=float, default=1.1)
    ap.add_argument("--min_samples", type=int, default=175)
    ap.add_argument("--max_samples", type=int, default=5000)
    ap.add_argument("--pval", type=float, default=0.05)
    ap.add_argument("--log2fc", type=float, default=0.5)
    ap.add_argument("--threshold_logic", type=str, default="and", choices=["and", "or"])
    ap.add_argument("--min_targets", type=int, default=25)
    ap.add_argument("--cap_edges", type=int, default=None)
    ap.add_argument("--verbosity", type=int, default=0)
    ap.add_argument("--display", action="store_true")
    ap.add_argument("--skip_mcmc", action="store_true", default=False)
    ap.add_argument("--entities", type=str, default=None, help="Path to .entities file for Entrez->symbol annotation")

    args = ap.parse_args()

    params = ORNORInferenceParams(
        seed=args.seed,
        chains=args.chains,
        gr_level=args.gr_level,
        min_samples=args.min_samples,
        max_samples=args.max_samples,
        p_value_threshold=args.pval,
        log2fc_threshold=args.log2fc,
        threshold_logic=args.threshold_logic,
        min_evidence_targets=args.min_targets,
        cap_edges=args.cap_edges,
        verbosity=args.verbosity,
        skip_mcmc=args.skip_mcmc,
    )

    print(f"[ORNOR runner] signature: {args.signature}")
    print(f"[ORNOR runner] network:   {args.network}")
    print(f"[ORNOR runner] out_edges: {args.out_edges}")
    print(f"[ORNOR runner] out_tfs:   {args.out_tfs}")
    print(f"[ORNOR runner] top_edges: {args.top_edges}")
    print(f"[ORNOR runner] params: {params}")

    t0 = time.time()
    edges_df, tf_df = run_ornor(args.signature, args.network, params, display=args.display)
    dt = time.time() - t0
    print(f"[ORNOR runner] finished in {dt:.2f} seconds")

    if tf_df is None or tf_df.empty:
        print("[ORNOR runner] WARNING: TF table is empty.")
        tf_df = pd.DataFrame(columns=["TF", "symbol", "X", "T"])
    else:
        # Annotate with gene symbols if entities file provided (or auto-detect)
        entities_path = args.entities
        if not entities_path:
            # Auto-detect alongside the network file, then fall back to resources/
            net_dir = Path(args.network).parent
            root_dir = Path(__file__).resolve().parents[3]
            candidates = [
                (net_dir, "all_tissues.fixed.entities"),
                (net_dir, "all_tissues.entities"),
                (net_dir, "all_tissues.entities.clean.tsv"),
                (root_dir / "resources", "tcCHIP_uid_map.tsv"),
            ]
            for d, name in candidates:
                fp = d / name
                if fp.exists():
                    entities_path = str(fp)
                    break
        if entities_path:
            sym_map = _load_entrez_symbol_map(entities_path)
            if sym_map:
                def _resolve_symbol(tf_id):
                    s = str(tf_id).strip()
                    name = sym_map.get(s, "")
                    # Treat name as unmapped if it equals the numeric ID
                    return name if name and name != s else ""
                tf_df["symbol"] = tf_df["TF"].astype(str).map(_resolve_symbol)
    tf_df.to_csv(args.out_tfs, sep="\t", index=False)

    if edges_df is None or edges_df.empty:
        print("[ORNOR runner] WARNING: No edges produced.")
        edges_out = pd.DataFrame(columns=["source", "target", "score"])
    else:
        edges_df = edges_df.copy()
        # Sort by X_tf descending (positive activators first) so the top-ranked
        # TFs' edges appear first in the capped file — required for the network
        # graph which filters edges by displayed TFs.
        if "X_tf" in edges_df.columns:
            edges_out = edges_df.sort_values("X_tf", ascending=False).head(int(args.top_edges))
        else:
            edges_df["abs_score"] = edges_df["score"].abs()
            edges_out = edges_df.sort_values("abs_score", ascending=False).head(int(args.top_edges))
        edges_out = edges_out[["source", "target", "score"]].copy()

    edges_out.to_csv(args.out_edges, index=False)

    print(f"[ORNOR runner] wrote edges: {args.out_edges} (rows={len(edges_out)})")
    print(f"[ORNOR runner] wrote tfs:   {args.out_tfs} (rows={len(tf_df)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
