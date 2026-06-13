"""
drugcentral.py
==============
Phase 2+4: Add DrugCentral structures and synonyms to the hub.

DrugCentral is licensed under CC BY-SA 4.0.
Source: https://drugcentral.org

SQL dump format (PostgreSQL COPY blocks, tab-delimited):
  structures  cols: cd_id(0), id(3), name(9), smiles(26), inchikey(29)
  synonyms    cols: syn_id(0), id(1), name(2), preferred_name(3)
  identifier  cols: id(0), identifier(1), id_type(2), struct_id(3)
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from hub import CompoundHub

_LICENSE = "cc-by-sa-4.0"
_SOURCE = "drugcentral"
_NULL = "\\N"

# Identifier types to import as xrefs
_XREF_TYPES = {
    "ChEMBL_ID": "CHEMBL_ID",
    "PUBCHEM_CID": "PUBCHEM_CID",
    "DRUGBANK_ID": "DRUGBANK_ID",
    "UNII": "UNII",
    "CHEBI": "CHEBI",
    "KEGG_DRUG": "KEGG_DRUG",
    "INN_ID": "INN_ID",
}


def _parse_sql_tables(filepath: str):
    """Stream the SQL dump, extract structures / synonyms / identifier tables.

    Returns three dicts keyed by DrugCentral compound id (string):
      structures  → {name, smiles, inchikey}
      synonyms    → [name, ...]          (preferred name first if present)
      identifiers → [(hub_id_type, value), ...]
    """
    structures: dict[str, dict] = {}
    synonyms_raw: dict[str, list] = {}   # dc_id -> [(name, is_preferred)]
    identifiers: dict[str, list] = {}    # dc_id -> [(hub_type, value)]

    current: str | None = None

    with open(filepath, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.rstrip("\n")

            if s.startswith("COPY public.structures "):
                current = "structures"
                continue
            if s.startswith("COPY public.synonyms "):
                current = "synonyms"
                continue
            if s.startswith("COPY public.identifier "):
                current = "identifier"
                continue
            if s == "\\.":
                current = None
                continue
            if s.startswith("COPY "):
                current = None
                continue

            if current is None:
                continue

            parts = s.split("\t")

            if current == "structures":
                if len(parts) < 30:
                    continue
                dc_id = parts[3]
                if not dc_id or dc_id == _NULL:
                    continue
                name = parts[9] if parts[9] != _NULL else ""
                smiles = parts[26] if parts[26] != _NULL else ""
                inchikey = parts[29] if parts[29] != _NULL else ""
                status = parts[30] if len(parts) > 30 and parts[30] != _NULL else ""
                structures[dc_id] = {"name": name, "smiles": smiles,
                                     "inchikey": inchikey, "status": status}

            elif current == "synonyms":
                if len(parts) < 3:
                    continue
                dc_id = parts[1]
                name = parts[2]
                preferred = parts[3] if len(parts) > 3 else _NULL
                if not dc_id or dc_id == _NULL or not name or name == _NULL:
                    continue
                synonyms_raw.setdefault(dc_id, []).append((name, preferred == "1"))

            elif current == "identifier":
                if len(parts) < 4:
                    continue
                id_type = parts[2]
                id_val = parts[1]
                dc_id = parts[3]
                hub_type = _XREF_TYPES.get(id_type)
                if not hub_type:
                    continue
                if not dc_id or dc_id == _NULL or not id_val or id_val == _NULL:
                    continue
                identifiers.setdefault(dc_id, []).append((hub_type, id_val))

    # Sort synonyms so preferred name comes first
    synonyms_out: dict[str, list[str]] = {}
    for dc_id, pairs in synonyms_raw.items():
        preferred = [n for n, p in pairs if p]
        other = [n for n, p in pairs if not p]
        synonyms_out[dc_id] = preferred + other

    return structures, synonyms_out, identifiers


def add_drugcentral(hub: CompoundHub, filepath: str) -> int:
    """Add DrugCentral structures, synonyms, and xrefs to the hub.

    Resolution strategy per compound:
      1. add_structure(smiles, inchikey) — deduplicates by InChIKey; returns
         existing node if already present, creates new node if novel.
      2. No structure → resolve by xrefs from identifier table.
      3. No xref match → resolve by normalized name.
      4. Still None → create orphan node.

    Returns the number of orphan nodes created.
    """
    if not os.path.exists(filepath):
        print(f"[Phase 2+4 / DrugCentral] not found, skipping: {filepath}")
        return 0

    print("  Parsing DrugCentral SQL dump …")
    structures, synonyms, identifiers = _parse_sql_tables(filepath)
    print(f"    structures={len(structures):,}  "
          f"synonym_groups={len(synonyms):,}  "
          f"identifier_groups={len(identifiers):,}")

    linked = new_struct = created_orphan = added_syns = 0

    for dc_id, struct in structures.items():
        names = synonyms.get(dc_id, [])
        pref_name = names[0] if names else struct["name"] or None
        smiles = struct["smiles"]
        inchikey = struct["inchikey"]
        status = struct.get("status", "")
        xrefs = identifiers.get(dc_id, [])

        sid: str | None = None

        if smiles or inchikey:
            before = len(hub._compounds)
            sid = hub.add_structure(
                smiles=smiles or None,
                inchikey=inchikey or None,
                source=_SOURCE,
                preferred_name=pref_name,
            )
            if sid is not None:
                if len(hub._compounds) > before:
                    new_struct += 1
                else:
                    linked += 1

        if sid is None and xrefs:
            sid = hub.resolve(xrefs=xrefs)
            if sid:
                linked += 1

        if sid is None and names:
            sid = hub.resolve(names=names[:5])
            if sid:
                linked += 1

        if sid is None:
            if not names and not struct["name"]:
                continue
            sid = hub.add_orphan(source=_SOURCE, preferred_name=pref_name)
            created_orphan += 1

        # Synonyms
        for name in names:
            if name and len(name) > 2:
                hub.add_synonym(sid, name, source=_SOURCE, license=_LICENSE)
                added_syns += 1

        # DrugCentral compound ID as xref
        hub.add_xref(sid, "DRUGCENTRAL", dc_id, source=_SOURCE, license=_LICENSE)

        # Approval status (ONP = Rx, OFP = Rx/OTC, OFM = OTC)
        if status:
            hub.add_xref(sid, "DC_STATUS", status, source=_SOURCE, license=_LICENSE)

        # Other identifier xrefs
        for hub_type, val in xrefs:
            hub.add_xref(sid, hub_type, val, source=_SOURCE, license=_LICENSE)

    print(f"[Phase 2+4 / DrugCentral]  linked={linked}  "
          f"new_structures={new_struct}  orphans={created_orphan}  "
          f"synonyms_added={added_syns}")

    return created_orphan
