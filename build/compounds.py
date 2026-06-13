"""
compounds.py
============
Phase 2: Build canonical compound nodes from structure-bearing sources.

Primary structure source: LINCS compoundinfo_beta (has canonical_smiles + inchi_key).
Additional sources can be added here as they become available (TTD SMILES download,
Drug Repurposing Hub samples, DrugCentral, B3DB, etc.).

Each loader returns the raw DataFrame for downstream use in xrefs.py / synonyms_build.py.
"""

from __future__ import annotations

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

import pandas as pd
from hub import CompoundHub

_SENTINEL = {"", "nan", '""', "none", "n/a", "na", "null"}


def _clean(val) -> str:
    s = str(val).strip().strip('"')
    return "" if s.lower() in _SENTINEL else s


def load_lincs(hub: CompoundHub, filepath: str) -> pd.DataFrame:
    """Load LINCS compoundinfo_beta.txt. Adds a canonical node for every record
    that has a canonical_smiles or inchi_key. Returns the raw DataFrame.

    Columns used: pert_id, cmap_name, canonical_smiles, inchi_key, compound_aliases
    """
    df = pd.read_csv(filepath, sep="\t", low_memory=False, dtype=str)
    df = df.fillna("")

    added = skipped = 0
    for _, row in df.iterrows():
        smiles = _clean(row.get("canonical_smiles", ""))
        inchikey = _clean(row.get("inchi_key", ""))
        name = _clean(row.get("cmap_name", "")) or _clean(row.get("pert_id", ""))

        if smiles:
            sid = hub.add_structure(smiles=smiles, source="lincs_2020",
                                    preferred_name=name or None)
        elif inchikey:
            sid = hub.add_structure(inchikey=inchikey, source="lincs_2020",
                                    preferred_name=name or None)
        else:
            skipped += 1
            continue

        if sid:
            added += 1

    print(f"[Phase 2 / LINCS]  nodes added={added}  skipped_no_structure={skipped}")
    return df


def load_repurposing_hub_samples(hub: CompoundHub, filepath: str) -> pd.DataFrame:
    """Optional: Load Drug Repurposing Hub samples file if available.
    Expected columns: broad_id, pert_iname, clinical_phase, moa, target,
                      disease_area, indication, smiles, InChIKey, pubchem_cid.
    """
    if not os.path.exists(filepath):
        print(f"[Phase 2 / RepHub] not found, skipping: {filepath}")
        return pd.DataFrame()

    df = pd.read_csv(filepath, dtype=str).fillna("")
    added = skipped = 0
    for _, row in df.iterrows():
        smiles = _clean(row.get("smiles", ""))
        inchikey = _clean(row.get("InChIKey", ""))
        name = _clean(row.get("pert_iname", ""))

        if smiles:
            sid = hub.add_structure(smiles=smiles, source="repurposing_hub_samples",
                                    preferred_name=name or None)
        elif inchikey:
            sid = hub.add_structure(inchikey=inchikey, source="repurposing_hub_samples",
                                    preferred_name=name or None)
        else:
            skipped += 1
            continue

        if sid:
            added += 1

    print(f"[Phase 2 / RepHub] nodes added={added}  skipped_no_structure={skipped}")
    return df
