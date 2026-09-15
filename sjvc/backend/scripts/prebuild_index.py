#!/usr/bin/env python3
"""Pre-build the SQLite index for a GENCODE release from a local GTF.

The path is also recorded in gencode_sources.json so later runs (even after the
cache is cleared) use the local file instead of downloading.

Usage:
    python scripts/prebuild_index.py human v29 ~/Work/SJ/Data/gencode.v29.annotation.gtf
    python scripts/prebuild_index.py human v29 <path> --force   # rebuild
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app.services.gencode import ensure_release, register_local_source  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("species", choices=["human", "mouse"])
    ap.add_argument("release", help="e.g. v29 (human) or vM25 (mouse)")
    ap.add_argument("gtf", type=Path, help="local .gtf or .gtf.gz")
    ap.add_argument("--force", action="store_true", help="rebuild even if cached")
    args = ap.parse_args()

    gtf = args.gtf.expanduser().resolve()
    if not gtf.exists():
        ap.error(f"no such file: {gtf}")

    register_local_source(args.species, args.release, gtf)
    print(f"registered {args.species}/{args.release} -> {gtf}")

    t0 = time.time()
    ann = ensure_release(args.species, args.release, force=args.force)
    dt = time.time() - t0
    print(f"index ready: {ann.db_path}  ({dt:.1f}s)")

    # quick sanity probe
    con = ann._con()
    try:
        n_genes = con.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
        n_feats = con.execute("SELECT COUNT(*) FROM feats").fetchone()[0]
    finally:
        con.close()
    print(f"  {n_genes:,} genes, {n_feats:,} exon/CDS rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
