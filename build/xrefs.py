"""
xrefs.py
========
Phase 3: Attach external cross-references to canonical compound nodes.

BRD is just one id_type here - no special status. Every xref row carries
source + license for the dual-build (Phase 6).
"""

from __future__ import annotations

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

import pandas as pd
from hub import CompoundHub
from normalize import normalize_name

_SENTINEL = {"", "nan", '""', "none", "n/a", "na", "null"}


def _clean(val) -> str:
    s = str(val).strip().strip('"')
    return "" if s.lower() in _SENTINEL else s


def add_lincs_xrefs(hub: CompoundHub, df_lincs: pd.DataFrame) -> None:
    """Attach BRD (pert_id) xrefs from the LINCS DataFrame.

    Resolves nodes by re-calling hub.add_structure() which is idempotent and
    returns the existing syndra_id when the InChIKey is already present.
    This avoids the chicken-and-egg of resolving by BRD before BRD is added.
    """
    added = 0
    for _, row in df_lincs.iterrows():
        smiles = _clean(row.get("canonical_smiles", ""))
        inchikey = _clean(row.get("inchi_key", ""))
        brd = _clean(row.get("pert_id", ""))

        # Resolve by structure (idempotent - returns existing node if InChIKey known)
        sid = hub.add_structure(
            smiles=smiles or None,
            inchikey=inchikey or None,
            source="lincs_2020",
        )
        if sid and brd:
            hub.add_xref(sid, "BRD", brd, source="lincs_2020",
                         license="custom_non_commercial")
            added += 1

    print(f"[Phase 3 / LINCS xrefs]  BRD entries linked={added}")


def add_prism_xrefs(hub: CompoundHub, filepath: str) -> int:
    """Link PubChem CIDs from PRISM to existing nodes by name or CID.

    Must run AFTER synonyms are populated (Phase 4) so name-based resolution works.
    """
    if not os.path.exists(filepath):
        print(f"[Phase 3 / PRISM xrefs] not found, skipping: {filepath}")
        return 0

    df = pd.read_csv(filepath, dtype=str).fillna("")
    linked = unlinked = 0

    for _, row in df.iterrows():
        drug_name = _clean(row.get("PRISM_drug_name", ""))
        pubchem = _clean(row.get("PubChem_CID", ""))

        sid = hub.resolve(
            xrefs=[("PUBCHEM_CID", pubchem)] if pubchem else None,
            names=[drug_name] if drug_name else None,
        )
        if sid is None:
            unlinked += 1
            continue

        if pubchem:
            hub.add_xref(sid, "PUBCHEM_CID", pubchem, source="prism",
                         license="cc-by-4.0")
        linked += 1

    print(f"[Phase 3 / PRISM xrefs]  linked={linked}  unlinked={unlinked}")
    return linked


def add_ttd_xrefs(hub: CompoundHub, df_ttd: pd.DataFrame) -> None:
    """Attach TTD IDs to nodes resolved by name.

    Must run AFTER synonyms are populated (Phase 4) so name-based resolution works.
    """
    linked = 0
    for _, row in df_ttd.iterrows():
        ttd_id = _clean(str(row.get("ttd_id", "")))
        drug_name = _clean(str(row.get("drug_name", "")))

        sid = hub.resolve(names=[drug_name] if drug_name else None)
        if sid and ttd_id:
            hub.add_xref(sid, "TTD", ttd_id, source="ttd", license="custom_academic")
            linked += 1

    print(f"[Phase 3 / TTD xrefs]  linked={linked}")
