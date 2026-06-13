"""
hub.py
======
The SynDRA canonical identifier hub (UniChem-style).

Three linked tables:
  compounds  - one canonical node per compound, keyed by `syndra_id`
  xrefs      - external IDs (BRD, PubChem CID, ChEMBL, UNII, TTD, DrugBank, ...)
  synonyms   - names, with the normalized match key

`syndra_id` is the canonical key for everything downstream. BRD is just one
`id_type` in xrefs. Nothing is dropped for lacking a BRD ID (orphan support).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

from normalize import normalize_name
from structure import standardize


@dataclass
class _Compound:
    syndra_id: str
    inchikey: Optional[str]
    inchikey_skeleton: Optional[str]
    standard_inchi: Optional[str]
    canonical_smiles: Optional[str]
    preferred_name: Optional[str]
    has_structure: bool
    first_source: Optional[str]


class CompoundHub:
    def __init__(self, id_prefix: str = "SYN", id_width: int = 7):
        self._id_prefix = id_prefix
        self._id_width = id_width
        self._counter = 0

        self._compounds: dict[str, _Compound] = {}
        self._inchikey_to_id: dict[str, str] = {}
        self._xref_to_id: dict[tuple[str, str], str] = {}
        self._synonym_to_ids: dict[str, set[str]] = {}

        self._xref_rows: dict[tuple, dict] = {}
        self._synonym_rows: dict[tuple, dict] = {}

        self.ambiguous_names: dict[str, set[str]] = {}

    def _mint(self) -> str:
        self._counter += 1
        return f"{self._id_prefix}{self._counter:0{self._id_width}d}"

    def add_structure(
        self,
        smiles: Optional[str] = None,
        inchikey: Optional[str] = None,
        source: Optional[str] = None,
        preferred_name: Optional[str] = None,
    ) -> Optional[str]:
        """Create or fetch the canonical node for a structure. Dedups by full
        standardized InChIKey. Returns syndra_id or None if unusable."""
        std = standardize(smiles) if smiles else None
        key = std.inchikey if std else (str(inchikey).strip() if inchikey else None)
        if not key:
            return None

        if key in self._inchikey_to_id:
            sid = self._inchikey_to_id[key]
            node = self._compounds[sid]
            if preferred_name and not node.preferred_name:
                node.preferred_name = preferred_name
            if std and not node.has_structure:
                node.inchikey = std.inchikey
                node.inchikey_skeleton = std.inchikey_skeleton
                node.standard_inchi = std.standard_inchi
                node.canonical_smiles = std.canonical_smiles
                node.has_structure = True
            return sid

        sid = self._mint()
        self._compounds[sid] = _Compound(
            syndra_id=sid,
            inchikey=key,
            inchikey_skeleton=(std.inchikey_skeleton if std else key.split("-")[0]),
            standard_inchi=(std.standard_inchi if std else None),
            canonical_smiles=(std.canonical_smiles if std else None),
            preferred_name=preferred_name,
            has_structure=bool(std),
            first_source=source,
        )
        self._inchikey_to_id[key] = sid
        return sid

    def add_orphan(self, source: Optional[str] = None,
                   preferred_name: Optional[str] = None) -> str:
        """Create a name-only node (no structure). Records that resolve to no
        structure AND no existing node are NOT dropped - this is the fix for
        the old BRD-only drop behavior."""
        sid = self._mint()
        self._compounds[sid] = _Compound(
            syndra_id=sid, inchikey=None, inchikey_skeleton=None,
            standard_inchi=None, canonical_smiles=None,
            preferred_name=preferred_name, has_structure=False, first_source=source,
        )
        return sid

    def add_xref(self, syndra_id: str, id_type: str, id_value: str,
                 source: Optional[str] = None, license: Optional[str] = None) -> None:
        id_type = str(id_type).strip().upper()
        id_value = str(id_value).strip()
        if not id_value or id_value in ("NAN", "") or syndra_id not in self._compounds:
            return
        self._xref_to_id.setdefault((id_type, id_value), syndra_id)
        self._xref_rows[(syndra_id, id_type, id_value, source)] = {
            "syndra_id": syndra_id, "id_type": id_type, "id_value": id_value,
            "source": source, "license": license,
        }

    def add_synonym(self, syndra_id: str, raw: str, source: Optional[str] = None,
                    license: Optional[str] = None, synonym_type: Optional[str] = None) -> None:
        norm = normalize_name(raw)
        if not norm or syndra_id not in self._compounds:
            return
        ids = self._synonym_to_ids.setdefault(norm, set())
        ids.add(syndra_id)
        if len(ids) > 1:
            self.ambiguous_names[norm] = set(ids)
        self._synonym_rows[(syndra_id, norm, source)] = {
            "syndra_id": syndra_id, "synonym_raw": str(raw).strip(),
            "synonym_norm": norm, "synonym_type": synonym_type,
            "source": source, "license": license,
        }

    def resolve(self, xrefs: Optional[list[tuple[str, str]]] = None,
                names: Optional[list[str]] = None) -> Optional[str]:
        """Find an existing node by xref or normalized name. xref match wins.
        Returns syndra_id or None."""
        if xrefs:
            for id_type, id_value in xrefs:
                sid = self._xref_to_id.get(
                    (str(id_type).strip().upper(), str(id_value).strip()))
                if sid:
                    return sid
        if names:
            for n in names:
                norm = normalize_name(n)
                hit = self._synonym_to_ids.get(norm)
                if hit:
                    if len(hit) > 1:
                        self.ambiguous_names[norm] = set(hit)
                    return sorted(hit)[0]
        return None

    def export(self) -> dict[str, pd.DataFrame]:
        compounds = pd.DataFrame([c.__dict__ for c in self._compounds.values()])
        xrefs = pd.DataFrame(list(self._xref_rows.values()))
        synonyms = pd.DataFrame(list(self._synonym_rows.values()))
        return {"compounds": compounds, "xrefs": xrefs, "synonyms": synonyms}

    def export_redistributable(self, prohibited_sources=None,
                               prohibited_licenses=None) -> dict[str, pd.DataFrame]:
        """Tables with rows from redistribution-prohibited sources removed.
        Compound nodes are kept; only offending attribution rows are stripped."""
        prohibited_sources = {str(s).lower() for s in (prohibited_sources or [])}
        prohibited_licenses = {str(l).lower() for l in (prohibited_licenses or [])}

        def ok(row) -> bool:
            s = str(row.get("source") or "").lower()
            l = str(row.get("license") or "").lower()
            return s not in prohibited_sources and l not in prohibited_licenses

        out = self.export()
        for tbl in ("xrefs", "synonyms"):
            df = out[tbl]
            if not df.empty:
                out[tbl] = df[df.apply(ok, axis=1)].reset_index(drop=True)
        return out

    def summary(self) -> dict:
        n_struct = sum(1 for c in self._compounds.values() if c.has_structure)
        return {
            "compounds": len(self._compounds),
            "with_structure": n_struct,
            "orphans_no_structure": len(self._compounds) - n_struct,
            "xref_rows": len(self._xref_rows),
            "synonym_rows": len(self._synonym_rows),
            "ambiguous_names": len(self.ambiguous_names),
        }
