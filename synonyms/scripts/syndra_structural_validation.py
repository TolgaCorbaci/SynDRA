#!/usr/bin/env python3
"""
SynDRA structural validation (OFFLINE).

Validates synonym -> BROAD_drug_ID merges against chemical structure, using the
canonical_smiles already present in LINCS compoundinfo_beta. Salt/solvate forms
of the same parent compound are treated as identical via salt stripping, so that
e.g. "imatinib" and "imatinib mesylate" resolve to the same drug.

Canonical drug key = first (connectivity) block of the InChIKey of the
standardised PARENT (largest organic fragment, neutralised). This is robust to:
  - salts / solvates / counterions  (stripped to parent)
  - charge / protonation state       (neutralised)
  - stereochemistry / isotopes       (dropped by using the connectivity block)

No internet required. For external query drug lists that lack structures, resolve
names -> SMILES first (see companion script), then reuse parent_inchikey() here.

Requires: pandas, rdkit
"""
import argparse
from functools import lru_cache
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import inchi

RDLogger.DisableLog("rdApp.*")  # silence parse warnings
_LFC = rdMolStandardize.LargestFragmentChooser()
_UNCHARGER = rdMolStandardize.Uncharger()


@lru_cache(maxsize=None)
def parent_inchikey(smiles, connectivity_only=True):
    """SMILES -> standardised-parent InChIKey (or its connectivity block). None on failure."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = _LFC.choose(mol)        # strip counterions/solvents
        mol = _UNCHARGER.uncharge(mol)  # neutralise
        key = inchi.MolToInchiKey(mol)
    except Exception:
        return None
    if not key:
        return None
    return key.split("-")[0] if connectivity_only else key


def run(compound_info, synonym_table):
    comp = pd.read_csv(compound_info, sep="\t", dtype=str)
    syn = pd.read_csv(synonym_table, dtype=str)
    syn["synonyms"] = syn["synonyms"].str.lower().str.strip()

    # BROAD_drug_ID -> parent connectivity key (structural ground truth)
    comp["parent_key"] = comp["canonical_smiles"].map(parent_inchikey)
    brd2key = dict(zip(comp["pert_id"], comp["parent_key"]))
    resolved = comp["parent_key"].notna().sum()
    print(f"Compounds in catalog: {len(comp):,} | structures resolved to a parent key: "
          f"{resolved:,} ({resolved/len(comp)*100:.1f}%)")

    # ---- A. Precision of the merge: classify synonym->multiple-BRD collisions ----
    def is_junk(s):
        s = str(s)
        return s.replace(".", "").replace("-", "").replace(" ", "").isdigit() or len(s) <= 2

    g = (syn.dropna(subset=["synonyms", "BROAD_drug_ID"])
            .groupby("synonyms")["BROAD_drug_ID"].apply(lambda x: sorted(set(x))))
    collisions = g[g.map(len) > 1]
    benign = true_error = junk = unresolved = 0
    error_examples = []
    for s, brds in collisions.items():
        if is_junk(s):
            junk += 1
            continue
        keys = {brd2key.get(b) for b in brds}
        keys.discard(None)
        if not keys:
            unresolved += 1
        elif len(keys) == 1:
            benign += 1  # same parent -> salt/stereo variants; synonym legitimately matches both
        else:
            true_error += 1
            if len(error_examples) < 12:
                error_examples.append((s, brds))
    total = len(collisions)
    print(f"\n[A] Synonym collisions (1 synonym -> >1 BROAD ID): {total}")
    print(f"    junk-string (numeric/<=2char): {junk}")
    print(f"    benign  (same parent structure, i.e. salt/form variants): {benign}")
    print(f"    TRUE merge errors (different parent structures):          {true_error}")
    print(f"    unresolved (no structure available):                      {unresolved}")
    structurally_checkable = benign + true_error
    if structurally_checkable:
        prec = benign / structurally_checkable * 100
        print(f"    -> merge precision on structurally-checkable collisions: {prec:.1f}%")
    for s, brds in error_examples:
        print(f"       ERROR  '{s}' -> {brds}  parents={[brd2key.get(b) for b in brds]}")

    # ---- B. Salt/form redundancy: BRD IDs that collapse to one parent compound ----
    cc = comp.dropna(subset=["parent_key"])
    by_parent = cc.groupby("parent_key")["pert_id"].nunique()
    redundant = by_parent[by_parent > 1]
    print(f"\n[B] Catalog redundancy by parent structure:")
    print(f"    unique parent compounds: {by_parent.size:,}")
    print(f"    parents represented by >1 BROAD ID (salt/form duplicates): {redundant.size:,}")
    print(f"    BROAD IDs involved in those duplicates: {int(redundant.sum()):,}")
    print(f"    largest single-parent BRD cluster: {int(redundant.max()) if redundant.size else 0} BROAD IDs")

    # ---- C. Junk synonym load (immediate cleaning target) ----
    s = syn["synonyms"].dropna()
    numeric_only = s.map(lambda x: str(x).replace(".", "").replace("-", "").replace(" ", "").isdigit()).sum()
    len2 = (s.str.len() <= 2).sum()
    print(f"\n[C] Junk synonyms acting as match keys: numeric-only={numeric_only:,}, <=2 chars={len2:,}")

    return syn, comp


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--compound-info", default="data/compoundinfo_beta.tsv")
    p.add_argument("--synonyms", default="data/merged_200K_drug_synonyms.csv")
    args = p.parse_args()
    run(args.compound_info, args.synonyms)
