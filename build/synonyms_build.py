"""
synonyms_build.py
=================
Phase 4 + 5: Attach synonyms and handle orphan records.

Phase 4 - Direct: synonyms attached to structure-bearing LINCS records.
Phase 4 - Indirect: synonym-only sources (TTD, PRISM) resolved into
          existing nodes via shared xref or normalized name.
Phase 5 - Orphans: unresolvable name-only records are KEPT as provisional nodes
          (not dropped), fixing the old BRD-centric drop behavior.
"""

from __future__ import annotations

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

import pandas as pd
from hub import CompoundHub
from normalize import normalize_name, split_synonyms

_SENTINEL = {"", "nan", '""', "none", "n/a", "na", "null"}


def _clean(val) -> str:
    s = str(val).strip().strip('"')
    return "" if s.lower() in _SENTINEL else s


def _is_junk(name: str) -> bool:
    """Filter out numeric-only and ≤2-char synonyms (per legacy pipeline)."""
    n = name.strip()
    if len(n) <= 2:
        return True
    if n.replace(".", "").replace("-", "").replace(" ", "").isdigit():
        return True
    return False


# ---------------------------------------------------------------------------
# TTD parser
# ---------------------------------------------------------------------------

def parse_ttd(filepath: str) -> pd.DataFrame:
    """Parse TTD P1-04-Drug_synonyms.txt into a flat DataFrame with columns:
    ttd_id, drug_name, synonyms (list).

    The file has a text preamble followed by records in the format:
      {TTD_ID}\\t{field}\\t{value}
    where field is one of: TTDDRUID, DRUGNAME, SYNONYMS
    """
    records: dict[str, dict] = {}
    preamble_done = False

    with open(filepath, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            ttd_id, field, value = parts[0].strip(), parts[1].strip(), parts[2].strip()
            if not ttd_id.startswith("D"):
                continue
            rec = records.setdefault(ttd_id, {"ttd_id": ttd_id, "drug_name": "", "synonyms": []})
            if field == "DRUGNAME":
                rec["drug_name"] = value
            elif field == "SYNONYMS":
                if value:
                    rec["synonyms"].append(value)

    rows = list(records.values())
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Phase 4: Direct synonyms from LINCS
# ---------------------------------------------------------------------------

def add_lincs_synonyms(hub: CompoundHub, df_lincs: pd.DataFrame) -> None:
    """Add cmap_name and pipe-separated compound_aliases from LINCS as synonyms.

    Resolves nodes by structure (idempotent hub.add_structure call), since BRD xrefs
    are added in the same Phase 3 step and may not resolve by name yet.
    """
    added = 0
    for _, row in df_lincs.iterrows():
        smiles = _clean(row.get("canonical_smiles", ""))
        inchikey = _clean(row.get("inchi_key", ""))
        name = _clean(row.get("cmap_name", ""))

        # Resolve by structure - idempotent, returns existing node
        sid = hub.add_structure(
            smiles=smiles or None,
            inchikey=inchikey or None,
            source="lincs_2020",
        )
        if sid is None:
            continue

        if name and not _is_junk(name):
            hub.add_synonym(sid, name, source="lincs_2020",
                            license="custom_non_commercial", synonym_type="INN")
            added += 1

        aliases_raw = _clean(row.get("compound_aliases", ""))
        for alias in split_synonyms(aliases_raw):
            if alias and not _is_junk(alias):
                hub.add_synonym(sid, alias, source="lincs_2020",
                                license="custom_non_commercial")
                added += 1

    print(f"[Phase 4 / LINCS synonyms]  synonym rows added≈{added}")


# ---------------------------------------------------------------------------
# Phase 4: Indirect - TTD
# ---------------------------------------------------------------------------

def add_ttd_synonyms(hub: CompoundHub, df_ttd: pd.DataFrame) -> int:
    """Add TTD drug names and synonyms. Link to existing nodes where possible;
    create orphan nodes for unresolvable TTD entries (Phase 5)."""
    linked = created_orphan = 0

    for _, row in df_ttd.iterrows():
        ttd_id = _clean(str(row.get("ttd_id", "")))
        drug_name = _clean(str(row.get("drug_name", "")))
        synonyms = list(row.get("synonyms", []))

        all_names = [n for n in ([drug_name] + synonyms) if n and not _is_junk(n)]
        if not all_names:
            continue

        sid = hub.resolve(
            xrefs=[("TTD", ttd_id)] if ttd_id else None,
            names=all_names,
        )

        if sid is None:
            # Phase 5: keep as orphan
            sid = hub.add_orphan(source="ttd",
                                 preferred_name=drug_name or all_names[0])
            created_orphan += 1
        else:
            linked += 1

        for name in all_names:
            hub.add_synonym(sid, name, source="ttd", license="custom_academic")

    print(f"[Phase 4+5 / TTD]  linked={linked}  orphans_created={created_orphan}")
    return created_orphan


# ---------------------------------------------------------------------------
# Phase 4: Indirect - PRISM
# ---------------------------------------------------------------------------

def add_prism_synonyms(hub: CompoundHub, filepath: str) -> int:
    """Add PRISM drug names and PubChem synonyms. Link via PubChem CID or name;
    create orphans for unresolvable entries."""
    if not os.path.exists(filepath):
        print(f"[Phase 4 / PRISM] not found, skipping: {filepath}")
        return 0

    df = pd.read_csv(filepath, dtype=str).fillna("")
    linked = created_orphan = 0

    for _, row in df.iterrows():
        drug_name = _clean(row.get("PRISM_drug_name", ""))
        pubchem = _clean(row.get("PubChem_CID", ""))
        raw_syns = _clean(row.get("PubChem_synonyms", ""))

        # PubChem_synonyms uses pipe as delimiter
        pub_syns = [s.strip() for s in raw_syns.split("|") if s.strip()] if raw_syns else []
        all_names = [n for n in ([drug_name] + pub_syns) if n and not _is_junk(n)]

        sid = hub.resolve(
            xrefs=[("PUBCHEM_CID", pubchem)] if pubchem else None,
            names=all_names[:5],  # limit name-match candidates to first 5
        )

        if sid is None:
            if not all_names:
                continue
            sid = hub.add_orphan(source="prism", preferred_name=drug_name or all_names[0])
            created_orphan += 1
        else:
            linked += 1

        for name in all_names[:20]:  # cap at 20 synonyms per PRISM entry
            hub.add_synonym(sid, name, source="prism", license="cc-by-4.0")

    print(f"[Phase 4+5 / PRISM]  linked={linked}  orphans_created={created_orphan}")
    return created_orphan
