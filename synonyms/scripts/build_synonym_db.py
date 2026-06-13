#!/usr/bin/env python3
"""
build_synonym_db.py — reproducible SynDRA build pipeline.

Integrates drug synonyms from three sources into a single
    synonym -> {BROAD_drug_ID, TTD_drug_ID, PubChem_CID}
mapping, then propagates identifiers across sources through shared synonyms and
keeps the BROAD-anchored rows (LINCS/L1000 focus).

This replaces the exploratory `create_synonym_database.ipynb`:
  * all paths are repository-relative (configurable via CLI), no absolute paths
  * parses the RAW source files committed under ../input (incl. the raw TTD
    P1-04 long-format file with its text preamble)
  * synonyms are normalised consistently (lower-cased, stripped)
  * no dead/commented exploratory cells

Sources (see DATA_SOURCES.md for versions/licenses):
  LINCS   compoundinfo_beta.txt      pert_id    -> BROAD_drug_ID
  TTD     P1-04-Drug_synonyms.txt    TTD_drug_ID
  PRISM   PRISM_drug_synonyms.csv    PubChem_CID

Usage:
    python build_synonym_db.py
    python build_synonym_db.py --input-dir ../input --output data/merged_200K_drug_synonyms.csv
"""
import argparse
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd


def _norm(series):
    """Lower-case, strip, and blank-out empty/sentinel synonym strings."""
    s = series.astype("string").str.lower().str.strip()
    return s.replace({"": pd.NA, "nan": pd.NA, "none": pd.NA, "missing": pd.NA})


def load_lincs(path):
    df = pd.read_csv(path, sep="\t", engine="c")
    df = df.rename(columns={"pert_id": "BROAD_drug_ID"})
    # Both the primary name (cmap_name) and the alias list are valid synonyms.
    frames = [df[["BROAD_drug_ID", col]].rename(columns={col: "synonyms"})
              for col in ("cmap_name", "compound_aliases")]
    out = pd.concat(frames, ignore_index=True)
    out["synonyms"] = _norm(out["synonyms"])
    return (out.dropna(subset=["synonyms"]).drop_duplicates(subset=["synonyms"]))


def load_ttd(path):
    """Parse the raw TTD P1-04 long-format file (skips the text preamble).

    Layout after the preamble:  <TTD_ID>\\t<FIELD>\\t<VALUE>
    with FIELD in {TTDDRUID, DRUGNAME, SYNONYMS}. DRUGNAME and SYNONYMS values
    are both treated as synonyms for the drug.
    """
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) == 3 and parts[1] in ("DRUGNAME", "SYNONYMS"):
                rows.append((parts[0], parts[2]))
    df = pd.DataFrame(rows, columns=["TTD_drug_ID", "synonyms"])
    df["synonyms"] = _norm(df["synonyms"])
    return (df.dropna(subset=["synonyms"]).drop_duplicates(subset=["synonyms"]))


def load_prism(path):
    df = pd.read_csv(path)
    df["synonyms"] = df["PubChem_synonyms"].apply(
        lambda s: list(dict.fromkeys(x.strip().lower() for x in str(s).split("|")))
    )
    out = df.explode("synonyms")[["PubChem_CID", "synonyms"]]
    out["synonyms"] = _norm(out["synonyms"])
    return (out.dropna(subset=["synonyms"]).drop_duplicates(subset=["synonyms"]))


def groupwise_fill(df, group_keys, fill_columns):
    """Within rows sharing an identifier, forward/backward-fill the other IDs."""
    df = df.copy()
    for col in fill_columns:
        for key in group_keys:
            if col == key:
                continue
            valid = df[key].notna()
            df.loc[valid, col] = (
                df.loc[valid]
                  .groupby(key, group_keys=False)[col]
                  .transform(lambda x: x.ffill().bfill() if x.notna().any() else x)
            )
    return df


def iterative_groupwise_fill(df, group_keys, fill_columns, max_iter=10):
    df = df.copy()
    for _ in range(max_iter):
        before = df[fill_columns].isna().sum().sum()
        df = groupwise_fill(df, group_keys, fill_columns)
        if df[fill_columns].isna().sum().sum() >= before:
            break
    return df


def build(input_dir, output_path):
    input_dir = Path(input_dir)
    lincs = load_lincs(input_dir / "compoundinfo_beta.txt")
    ttd = load_ttd(input_dir / "P1-04-Drug_synonyms.txt")
    prism = load_prism(input_dir / "PRISM_drug_synonyms.csv")
    print(f"Loaded synonyms  LINCS={len(lincs):,}  "
          f"TTD={len(ttd):,}  PRISM={len(prism):,}")

    # Anchor source (BROAD-keyed), then 3-way outer merge on synonym
    broad = lincs.drop_duplicates(subset=["synonyms"])
    merged = reduce(lambda l, r: pd.merge(l, r, on="synonyms", how="outer"),
                    [broad, ttd, prism])
    merged["PubChem_CID"] = pd.to_numeric(merged["PubChem_CID"], errors="coerce")

    # Propagate identifiers across shared synonyms
    keys = ["PubChem_CID", "TTD_drug_ID", "BROAD_drug_ID"]
    merged = iterative_groupwise_fill(merged, keys, keys)

    # Keep BROAD-anchored rows (LINCS focus), dedupe, order columns
    merged = (merged.drop_duplicates()
                    .dropna(subset=["synonyms", "BROAD_drug_ID"])
                    .reset_index(drop=True))
    merged = merged[["BROAD_drug_ID", "synonyms", "TTD_drug_ID", "PubChem_CID"]]

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, index=False)
    print(f"\nWrote {len(merged):,} rows -> {output_path}")
    print(f"  unique synonyms : {merged['synonyms'].nunique():,}")
    print(f"  unique BROAD IDs: {merged['BROAD_drug_ID'].nunique():,}")
    print(f"  unique TTD IDs  : {merged['TTD_drug_ID'].nunique():,}")
    print(f"  unique PubChem  : {merged['PubChem_CID'].nunique():,}")
    return merged


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", default="../input",
                   help="folder with the four source files (default: ../input)")
    p.add_argument("--output", default="data/merged_200K_drug_synonyms.csv",
                   help="output CSV path (default: data/merged_200K_drug_synonyms.csv)")
    args = p.parse_args()
    build(args.input_dir, args.output)
