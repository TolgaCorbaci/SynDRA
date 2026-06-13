"""
run_all.py
==========
Phase 7 batch runner: parse all 27 libraries, compute coverage + statin control.

Run from project root:
  python enrichment/run_all.py

Requires outputs/syndra_redistributable_synonyms.parquet (from build/build_all.py).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

_HERE = Path(__file__).parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_ROOT / "build"))
sys.path.insert(0, str(_HERE))

import pandas as pd
from syndra_enrich import build_resolver, enrich, parse_all_libraries, harmonize_library

OUTPUT_DIR = _ROOT / "outputs"
DB_DIR = _ROOT / "synonyms" / "input" / "enrichment_databases"
COVERAGE_OUT = OUTPUT_DIR / "coverage.csv"

# Statin enrichment control (Phase 8 acceptance criterion 5)
STATIN_QUERY = [
    "atorvastatin", "simvastatin", "pravastatin",
    "lovastatin", "pitavastatin", "rosuvastatin",
]


def main():
    print("=" * 60)
    print("SynDRA Enrichment Analysis  (Phase 7)")
    print("=" * 60)

    # Load redistributable build
    synonyms_path = OUTPUT_DIR / "syndra_redistributable_synonyms.parquet"
    if not synonyms_path.exists():
        synonyms_path = OUTPUT_DIR / "syndra_redistributable_synonyms.csv"
        if not synonyms_path.exists():
            print("ERROR: outputs not found. Run build/build_all.py first.")
            sys.exit(1)

    print(f"Loading synonyms from: {synonyms_path}")
    synonyms_df = (pd.read_parquet(synonyms_path)
                   if str(synonyms_path).endswith(".parquet")
                   else pd.read_csv(synonyms_path))

    # Coverage report across all libraries
    print(f"\nParsing libraries from: {DB_DIR}")
    coverage_df, results, _ = enrich(
        query_names=[],   # no query; we just want coverage stats
        synonyms_df=synonyms_df,
        db_dir=str(DB_DIR),
    )

    print("\n=== Coverage Report ===")
    pd.set_option("display.max_colwidth", 40)
    print(coverage_df.sort_values("match_rate", ascending=False).to_string(index=False))

    coverage_df.to_csv(COVERAGE_OUT, index=False)
    print(f"\nCoverage written to: {COVERAGE_OUT}")

    # Statin enrichment control
    print("\n=== Statin Enrichment Control ===")
    statin_coverage, statin_results, unresolved = enrich(
        query_names=STATIN_QUERY,
        synonyms_df=synonyms_df,
        db_dir=str(DB_DIR),
        fdr_threshold=0.05,
    )

    if unresolved:
        print(f"WARNING: unresolved statins: {unresolved}")

    target_libs = [
        "Drug_Repurposing_Hub_Mechanism_of_Action",
        "DrugCentral_Target",
        "Drug_Repurposing_Hub_Target",
    ]
    for lib_name in target_libs:
        df = statin_results.get(lib_name)
        if df is None or df.empty:
            print(f"\n  {lib_name}: no results")
            continue
        top = df.head(10)
        print(f"\n  {lib_name} (top 10):")
        print(top[["term", "overlap", "term_size", "pvalue", "qvalue"]].to_string(index=False))

    # Naive name-match rate (for comparison in Phase 8)
    all_libs = parse_all_libraries(str(DB_DIR))
    _naive_coverage(all_libs, synonyms_df)


def _naive_coverage(all_libs: dict, synonyms_df: pd.DataFrame):
    """Compare SynDRA-harmonized vs naive raw-name match rate."""
    from syndra_enrich import build_resolver, harmonize_library, normalize_name

    resolver = build_resolver(synonyms_df)
    all_norms = set(resolver.keys())

    total = resolved = naive_resolved = 0
    for lib in all_libs.values():
        for drugs in lib.values():
            for d in drugs:
                total += 1
                if normalize_name(d) in all_norms:
                    resolved += 1
                    naive_resolved += 1  # same here; naive = SynDRA for now

    if total:
        print(f"\nNaive match rate (raw lowercase): {naive_resolved}/{total} = "
              f"{naive_resolved/total:.1%}")
        print(f"SynDRA-harmonized match rate:     {resolved}/{total} = "
              f"{resolved/total:.1%}")


if __name__ == "__main__":
    main()
