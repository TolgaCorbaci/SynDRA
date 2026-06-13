"""
recall.py
=========
Phase 8: Recall lift benchmark.

Compares fraction of 527-drug benchmark library resolved under:
  (a) Naive string match (LINCS cmap_name field only)
  (b) Structure-anchored SynDRA (syndra_redistributable synonyms)

Run from project root:
  python benchmark/recall.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_ROOT / "build"))
sys.path.insert(0, str(_ROOT / "enrichment"))

import pandas as pd
from normalize import normalize_name

INPUT = _ROOT / "synonyms" / "input"
OUTPUT = _ROOT / "outputs"

BENCHMARK_PATH = INPUT / "File_2_Drugname_library_527D.txt"
LINCS_PATH = INPUT / "compoundinfo_beta.txt"
NEW_SYNONYMS_PATH = OUTPUT / "syndra_redistributable_synonyms.parquet"


def load_benchmark(filepath: str) -> pd.DataFrame:
    """Load benchmark file (TSV: drug, InChiKey, synonyms, drugclass).
    Returns DataFrame with one row per benchmark drug.
    """
    return pd.read_csv(filepath, sep="\t", dtype=str).fillna("")


def _all_names_for_drug(row: pd.Series) -> list[str]:
    """Collect primary name + semicolon-separated synonyms for one benchmark drug."""
    names = [str(row.get("drug", "")).strip()]
    raw_syns = str(row.get("synonyms", "")).strip()
    if raw_syns:
        names += [s.strip() for s in raw_syns.split(";") if s.strip()]
    return [n for n in names if n]


def _matched_count(df_bench: pd.DataFrame, norm_set: set[str]) -> tuple[int, list[str]]:
    """Count benchmark drugs resolved into norm_set (any name or synonym)."""
    matched = unmatched = 0
    unmatched_drugs = []
    for _, row in df_bench.iterrows():
        names = _all_names_for_drug(row)
        if any(normalize_name(n) in norm_set for n in names):
            matched += 1
        else:
            unmatched += 1
            unmatched_drugs.append(row.get("drug", ""))
    return matched, unmatched_drugs


def recall_naive(df_bench: pd.DataFrame, lincs_path: str) -> tuple[int, list[str]]:
    """Method (a): naive match against LINCS cmap_name field only."""
    df = pd.read_csv(lincs_path, sep="\t", dtype=str).fillna("")
    lincs_norms = {normalize_name(n) for n in df["cmap_name"].tolist() if n}
    return _matched_count(df_bench, lincs_norms)


def recall_new_syndra(df_bench: pd.DataFrame, synonyms_path: str) -> tuple[int, list[str]]:
    """Method (c): new structure-anchored SynDRA (syndra_id-keyed synonyms)."""
    path = Path(synonyms_path)
    csv_path = path.with_suffix(".csv")

    if path.exists():
        df = pd.read_parquet(path)
    elif csv_path.exists():
        df = pd.read_csv(csv_path, dtype=str).fillna("")
    else:
        print("  [new SynDRA] outputs not found. Run build/build_all.py first.")
        return 0, list(df_bench["drug"])

    new_norms = {str(n) for n in df["synonym_norm"].tolist() if n}
    return _matched_count(df_bench, new_norms)


def main():
    print("=" * 60)
    print("SynDRA Recall Benchmark  (Phase 8)")
    print("=" * 60)

    df_bench = load_benchmark(str(BENCHMARK_PATH))
    N = len(df_bench)
    print(f"Benchmark set: {N} drugs")

    n_naive, miss_naive = recall_naive(df_bench, str(LINCS_PATH))
    n_new, miss_new = recall_new_syndra(df_bench, str(NEW_SYNONYMS_PATH))

    print("\n=== Recall Results ===")
    print(f"  (a) Naive (LINCS cmap_name):      {n_naive}/{N} = {n_naive/N:.1%}")
    print(f"  (b) SynDRA:                        {n_new}/{N} = {n_new/N:.1%}")
    print(f"\n  Lift (b) vs (a): {n_new - n_naive:+d} drugs ({(n_new-n_naive)/N:+.1%})")

    if miss_new:
        print(f"\n  Still unresolved ({len(miss_new)}):")
        for name in sorted(miss_new)[:20]:
            print(f"    {name}")
        if len(miss_new) > 20:
            print(f"    ... and {len(miss_new) - 20} more")


if __name__ == "__main__":
    main()
