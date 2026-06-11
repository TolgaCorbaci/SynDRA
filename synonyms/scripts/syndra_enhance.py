#!/usr/bin/env python3
"""
SynDRA enhancement (OFFLINE): upgrade the synonym resource from string-only to
structure-aware, and emit a structural de-duplication map for the L1000 catalog.

Outputs
  1. merged_synonyms_parent_augmented.csv
       original columns + `parent_inchikey` (salt-stripped connectivity key for
       the BROAD compound). Junk synonyms (numeric-only, <=2 chars) removed;
       structurally-ambiguous synonyms flagged in `ambiguous_flag`.
  2. brd_parent_consolidation.csv
       parent_inchikey -> the set of BROAD IDs that are salt/form/stereo variants
       of one parent compound. Lets a single query hit ALL replicate signatures.

Reuses parent_inchikey() from syndra_structural_validation.py.
Requires: pandas, rdkit
"""
import pandas as pd
from syndra_structural_validation import parent_inchikey

COMPOUND_INFO = "data/compoundinfo_beta.tsv"
SYNONYMS = "data/merged_200K_drug_synonyms.csv"
OUTDIR = "/mnt/user-data/outputs"


def is_junk(s):
    s = str(s)
    return s.replace(".", "").replace("-", "").replace(" ", "").isdigit() or len(s) <= 2


def main():
    comp = pd.read_csv(COMPOUND_INFO, sep="\t", dtype=str)
    syn = pd.read_csv(SYNONYMS, dtype=str)
    syn["synonyms"] = syn["synonyms"].str.lower().str.strip()

    comp["parent_inchikey"] = comp["canonical_smiles"].map(parent_inchikey)
    brd2key = dict(zip(comp["pert_id"], comp["parent_inchikey"]))

    # ---- 1. augment + clean the synonym table ----
    before = len(syn)
    syn = syn[~syn["synonyms"].map(is_junk)].copy()
    syn["parent_inchikey"] = syn["BROAD_drug_ID"].map(brd2key)

    # flag synonyms that resolve to >1 distinct parent structure (true ambiguity)
    syn_parents = (syn.dropna(subset=["synonyms", "parent_inchikey"])
                      .groupby("synonyms")["parent_inchikey"].nunique())
    ambiguous = set(syn_parents[syn_parents > 1].index)
    syn["ambiguous_flag"] = syn["synonyms"].isin(ambiguous)

    out1 = f"{OUTDIR}/merged_synonyms_parent_augmented.csv"
    syn.to_csv(out1, index=False)
    print(f"[1] {out1}")
    print(f"    rows {before:,} -> {len(syn):,} (junk synonyms removed: {before-len(syn):,})")
    print(f"    rows with a parent structure: {syn['parent_inchikey'].notna().sum():,} "
          f"({syn['parent_inchikey'].notna().mean()*100:.1f}%)")
    print(f"    synonyms flagged structurally ambiguous: {len(ambiguous)}")

    # ---- 2. structural consolidation map ----
    cc = comp.dropna(subset=["parent_inchikey"])
    cons = (cc.groupby("parent_inchikey")["pert_id"]
              .agg(n_brd_ids="nunique", brd_ids=lambda x: "; ".join(sorted(set(x)))))
    cons = cons.sort_values("n_brd_ids", ascending=False).reset_index()
    out2 = f"{OUTDIR}/brd_parent_consolidation.csv"
    cons.to_csv(out2, index=False)
    dup = cons[cons["n_brd_ids"] > 1]
    print(f"\n[2] {out2}")
    print(f"    unique parent compounds: {len(cons):,}")
    print(f"    parents with >1 BROAD ID: {len(dup):,}  "
          f"covering {int(dup['n_brd_ids'].sum()):,} BROAD IDs")
    print(f"    => {int(dup['n_brd_ids'].sum())/cc['pert_id'].nunique()*100:.0f}% of catalog "
          f"entries are salt/form/stereo variants of a smaller parent set")


if __name__ == "__main__":
    main()
