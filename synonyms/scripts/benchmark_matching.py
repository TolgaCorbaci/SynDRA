#!/usr/bin/env python3
"""
benchmark_matching.py — reproducible SynDRA benchmark.

Measures how many drugs in a query library can be matched to a BROAD_drug_ID
using (a) the full SynDRA synonym map vs (b) the LINCS `cmap_name` field alone
(the naive baseline). The gain is SynDRA's added recall.

Replaces the exploratory `synonym_matching.ipynb` with repository-relative paths
and no dead cells.

Query library format (tab-separated): columns `drug`, `synonyms` (`;`-separated),
optionally `InChiKey`, `drugclass`.

Usage:
    python benchmark_matching.py
    python benchmark_matching.py --library ../input/File_2_Drugname_library_527D.txt \
        --synonyms data/merged_200K_drug_synonyms.csv \
        --lincs ../input/compoundinfo_beta.txt
"""
import argparse
from pathlib import Path

import pandas as pd


def load_query_library(path):
    """Explode `drug` + `synonyms` into one lower-cased synonym per row."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    total_drugs = df["drug"].str.lower().str.strip().nunique()
    combined = (df["drug"].fillna("") + "; " + df["synonyms"].fillna("")).str.split(";")
    out = df.assign(synonyms=combined).explode("synonyms")
    out["synonyms"] = out["synonyms"].str.lower().str.strip()
    out["drug"] = out["drug"].str.lower().str.strip()
    out = (out[out["synonyms"].notna() & (out["synonyms"] != "")]
           .drop_duplicates(subset=["synonyms"]).reset_index(drop=True))
    return out, total_drugs


def match_rate(query, total_drugs, mapping, right_col):
    """Inner-join query synonyms to a mapping column; return (#drugs matched, %)."""
    mapping = mapping.copy()
    mapping[right_col] = mapping[right_col].str.lower().str.strip()
    merged = pd.merge(query, mapping, left_on="synonyms", right_on=right_col, how="inner")
    matched = merged["drug"].nunique()
    return matched, matched / total_drugs * 100 if total_drugs else 0.0


def run(library, synonyms_path, lincs_path):
    query, total = load_query_library(library)
    print(f"Query library: {total} drugs ({len(query):,} unique synonyms after explosion)\n")

    syndra = pd.read_csv(synonyms_path, dtype=str)
    s_matched, s_pct = match_rate(query, total, syndra, "synonyms")

    lincs = pd.read_csv(lincs_path, sep="\t", dtype=str)
    l_matched, l_pct = match_rate(query, total, lincs, "cmap_name")

    print(f"{'Method':<34}{'Drugs matched':>15}{'Match rate':>13}")
    print(f"{'-'*62}")
    print(f"{'LINCS cmap_name only (baseline)':<34}{l_matched:>15}{l_pct:>12.1f}%")
    print(f"{'SynDRA synonym map':<34}{s_matched:>15}{s_pct:>12.1f}%")
    print(f"{'-'*62}")
    print(f"{'Improvement':<34}{s_matched - l_matched:>15}{s_pct - l_pct:>+12.1f}%")
    return {"total": total, "baseline_pct": l_pct, "syndra_pct": s_pct,
            "improvement_pp": s_pct - l_pct}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--library", default="../input/File_2_Drugname_library_527D.txt")
    p.add_argument("--synonyms", default="data/merged_200K_drug_synonyms.csv")
    p.add_argument("--lincs", default="../input/compoundinfo_beta.txt")
    args = p.parse_args()
    run(args.library, args.synonyms, args.lincs)
