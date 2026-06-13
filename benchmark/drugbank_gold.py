"""
drugbank_gold.py
================
Phase 8: DrugBank gold-standard name -> InChIKey evaluation.

Uses the DrugBank 2021 academic file as a PRIVATE benchmark only.
DrugBank data is NEVER included in any redistributable output.

For each DrugBank drug name, checks whether SynDRA maps it to the same InChIKey.
Reports precision/recall.

Run from project root:
  python benchmark/drugbank_gold.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

_HERE = Path(__file__).parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_ROOT / "build"))

import pandas as pd
from normalize import normalize_name

DRUGBANK_VOCAB = _ROOT / "synonyms" / "input" / "Drugbank" / "drugbank vocabulary.csv"
DRUGBANK_TSV   = _ROOT / "synonyms" / "input" / "Drugbank" / "drugbank_2021.tsv"
OUTPUT = _ROOT / "outputs"


def load_drugbank_gold(vocab_path: str) -> pd.DataFrame:
    """Load DrugBank vocabulary CSV. Expected columns include:
    Common name, Standard InChI Key, Synonyms (comma-separated).
    """
    df = pd.read_csv(vocab_path, dtype=str).fillna("")
    return df


def load_syndra_compounds(output_dir: str) -> pd.DataFrame:
    """Load syndra_full_compounds (has inchikey column)."""
    for fname in ["syndra_full_compounds.parquet", "syndra_full_compounds.csv"]:
        path = Path(output_dir) / fname
        if path.exists():
            return (pd.read_parquet(path) if fname.endswith(".parquet")
                    else pd.read_csv(path, dtype=str).fillna(""))
    return pd.DataFrame()


def load_syndra_synonyms(output_dir: str) -> pd.DataFrame:
    """Load syndra_full_synonyms for name->syndra_id resolution."""
    for fname in ["syndra_full_synonyms.parquet", "syndra_full_synonyms.csv"]:
        path = Path(output_dir) / fname
        if path.exists():
            return (pd.read_parquet(path) if fname.endswith(".parquet")
                    else pd.read_csv(path, dtype=str).fillna(""))
    return pd.DataFrame()


def main():
    print("=" * 60)
    print("SynDRA DrugBank Gold-Standard Evaluation  (Phase 8)")
    print("  [DrugBank used under academic license for evaluation only]")
    print("  [DrugBank data NOT redistributed]")
    print("=" * 60)

    if not DRUGBANK_VOCAB.exists():
        print(f"DrugBank vocabulary not found: {DRUGBANK_VOCAB}")
        print("This evaluation requires the DrugBank academic license files.")
        sys.exit(1)

    compounds_df = load_syndra_compounds(str(OUTPUT))
    synonyms_df  = load_syndra_synonyms(str(OUTPUT))

    if compounds_df.empty or synonyms_df.empty:
        print("SynDRA outputs not found. Run build/build_all.py first.")
        sys.exit(1)

    # Build resolver: synonym_norm -> syndra_id -> inchikey
    syn_to_sid = {}
    for _, row in synonyms_df.iterrows():
        norm = str(row.get("synonym_norm", "")).strip()
        sid = str(row.get("syndra_id", "")).strip()
        if norm and sid:
            syn_to_sid[norm] = sid

    _SENTINELS = {"", "none", "nan", "null", "n/a"}
    sid_to_ik = {}
    for _, row in compounds_df.iterrows():
        sid = str(row.get("syndra_id", "")).strip()
        ik = str(row.get("inchikey", "")).strip()
        if sid and ik.lower() not in _SENTINELS:
            sid_to_ik[sid] = ik

    # Evaluate against DrugBank
    db_df = load_drugbank_gold(str(DRUGBANK_VOCAB))

    # Detect InChIKey column
    ik_col = next((c for c in db_df.columns if "inchi" in c.lower() and "key" in c.lower()),
                  None)
    name_col = next((c for c in db_df.columns if "name" in c.lower()), None)

    if not ik_col or not name_col:
        print(f"Unexpected DrugBank columns: {list(db_df.columns)}")
        sys.exit(1)

    tp = fp = fn = 0
    errors = []

    for _, row in db_df.iterrows():
        db_name = str(row[name_col]).strip()
        db_ik = str(row[ik_col]).strip()
        if not db_name or not db_ik:
            continue

        norm = normalize_name(db_name)
        sid = syn_to_sid.get(norm)

        if sid is None:
            fn += 1  # SynDRA didn't find it at all
            continue

        syndra_ik = sid_to_ik.get(sid, "")
        if not syndra_ik:
            fn += 1
            continue

        if syndra_ik == db_ik:
            tp += 1
        else:
            fp += 1
            errors.append((db_name, db_ik, syndra_ik))

    total = tp + fp + fn
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    print(f"\n=== DrugBank Gold Standard Results ===")
    print(f"  DrugBank entries evaluated: {total}")
    print(f"  True positive  (correct InChIKey): {tp}")
    print(f"  False positive (wrong InChIKey):   {fp}")
    print(f"  False negative (not found):         {fn}")
    print(f"\n  Precision: {precision:.3f}")
    print(f"  Recall:    {recall:.3f}")
    print(f"  F1:        {f1:.3f}")

    if errors:
        print(f"\n  First 10 InChIKey mismatches:")
        for name, db_ik, syndra_ik in errors[:10]:
            print(f"    {name!r:30s}  DrugBank={db_ik}  SynDRA={syndra_ik}")

    print("\nNote: DrugBank used under academic license for evaluation; not redistributed.")


if __name__ == "__main__":
    main()
