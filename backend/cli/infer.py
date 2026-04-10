from __future__ import annotations

import argparse
from pathlib import Path

from backend.app.schema import SignatureInput, NetworkInput, RunConfig
from backend.app.service import InferenceService
from backend.engines.cie.adapter import CIEEngine
from backend.engines.ornor.adapter import ORNOREngine


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=["CIE", "ORNOR"], required=True)
    ap.add_argument("--signature", required=True)
    ap.add_argument("--network", required=True)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--results_root", default="results")
    args = ap.parse_args()

    service = InferenceService(results_root=Path(args.results_root))

    sig = SignatureInput(path=Path(args.signature))
    net = NetworkInput(path=Path(args.network))

    if args.method == "CIE":
        engine = CIEEngine()
    else:
        engine = ORNOREngine(python_bin="python")

    out = service.run(
        engine=engine,
        signature=sig,
        network=net,
        config=RunConfig(method=args.method, seed=args.seed, params={}),
    )

    print("DONE")
    print("edges:", out.edges_csv)
    print("regs :", out.regulators_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
