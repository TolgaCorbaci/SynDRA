"""
hub.py
======
The SynDRA canonical identifier hub (UniChem-style).

Replaces the flat `synonym <-> BRD` table with three linked tables:
  compounds  - one canonical node per compound, keyed by `syndra_id`,
               backed by a standardized InChIKey (from structure.py)
  xrefs      - external IDs (BRD, PubChem CID, ChEMBL, UNII, TTD, DrugBank, ...)
  synonyms   - names, with the normalized match key (from normalize.py)

`syndra_id` is the canonical key for everything downstream. BRD is just one
`id_type` in xrefs - nothing is dropped for lacking it (orphan support below).

Intended build order (see SYNDRA_IMPLEMENTATION_GUIDE.md):
  Phase 2: add_structure(...) for every structure-bearing record   -> nodes
  Phase 3: add_xref(...)                                            -> xrefs
  Phase 4: add_synonym(...) directly, or resolve(...) then attach   -> synonyms
  Phase 5: add_orphan(...) for name-only records                   -> kept, flagged
  Phase 6: export() / export_redistributable(...)                  -> outputs

The real source-specific loaders (Repurposing Hub samples, TTD, merged_data,
B3DB, DrugCentral, KatDB, PRISM synonyms, ...) are left to the agent; they just
call these methods. The hard, error-prone parts (dedup, linking, orphan keep,
license split) live here so they're consistent across sources.
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
        # dedup / linking indexes
        self._inchikey_to_id: dict[str, str] = {}
        self._xref_to_id: dict[tuple[str, str], str] = {}    # (id_type, id_value) -> syndra_id
        self._synonym_to_ids: dict[str, set[str]] = {}       # synonym_norm -> {syndra_id}

        # row stores (keyed so re-adding the same fact is idempotent)
        self._xref_rows: dict[tuple, dict] = {}
        self._synonym_rows: dict[tuple, dict] = {}

        # names that resolved to more than one node during the build (audit these)
        self.ambiguous_names: dict[str, set[str]] = {}

    # ---- id minting -------------------------------------------------
    def _mint(self) -> str:
        self._counter += 1
        return f"{self._id_prefix}{self._counter:0{self._id_width}d}"

    # ---- nodes ------------------------------------------------------
    def add_structure(
        self,
        smiles: Optional[str] = None,
        inchikey: Optional[str] = None,
        source: Optional[str] = None,
        preferred_name: Optional[str] = None,
    ) -> Optional[str]:
        """Create or fetch the canonical node for a structure. Dedups by full
        standardized InChIKey. Provide `smiles` (preferred; it gets standardized)
        or a precomputed `inchikey`. Returns the syndra_id, or None if unusable.

        If a structure arrives for a name-only orphan that already exists under
        the same InChIKey, the orphan is upgraded in place.
        """
        std = standardize(smiles) if smiles else None
        key = std.inchikey if std else (str(inchikey).strip() if inchikey else None)
        if not key:
            return None

        if key in self._inchikey_to_id:
            sid = self._inchikey_to_id[key]
            node = self._compounds[sid]
            if preferred_name and not node.preferred_name:
                node.preferred_name = preferred_name
            if std and not node.has_structure:          # upgrade orphan
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
        """Create a name-only node (no structure). Use for records that resolve
        to no structure AND no existing node - so they are NOT dropped (this is
        the explicit fix for the old BRD-only drop behavior)."""
        sid = self._mint()
        self._compounds[sid] = _Compound(
            syndra_id=sid, inchikey=None, inchikey_skeleton=None,
            standard_inchi=None, canonical_smiles=None,
            preferred_name=preferred_name, has_structure=False, first_source=source,
        )
        return sid

    # ---- xrefs ------------------------------------------------------
    def add_xref(self, syndra_id: str, id_type: str, id_value: str,
                 source: Optional[str] = None, license: Optional[str] = None) -> None:
        id_type = str(id_type).strip().upper()
        id_value = str(id_value).strip()
        if not id_value or syndra_id not in self._compounds:
            return
        self._xref_to_id.setdefault((id_type, id_value), syndra_id)
        self._xref_rows[(syndra_id, id_type, id_value, source)] = {
            "syndra_id": syndra_id, "id_type": id_type, "id_value": id_value,
            "source": source, "license": license,
        }

    # ---- synonyms ---------------------------------------------------
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

    # ---- indirect linking (Phase 4) --------------------------------
    def resolve(self, xrefs: Optional[list[tuple[str, str]]] = None,
                names: Optional[list[str]] = None) -> Optional[str]:
        """Find an existing node that shares any given xref (id_type, id_value)
        or any normalized name. Returns a syndra_id or None. Use to attach
        synonym-only sources INTO existing nodes before creating an orphan.

        xref match wins over name match (structurally grounded). If a name is
        ambiguous (maps to >1 node), returns one deterministically and records
        it in `ambiguous_names` for auditing.
        """
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

    # ---- export -----------------------------------------------------
    def export(self) -> dict[str, pd.DataFrame]:
        compounds = pd.DataFrame([c.__dict__ for c in self._compounds.values()])
        xrefs = pd.DataFrame(list(self._xref_rows.values()))
        synonyms = pd.DataFrame(list(self._synonym_rows.values()))
        return {"compounds": compounds, "xrefs": xrefs, "synonyms": synonyms}

    def export_redistributable(self, prohibited_sources=None,
                               prohibited_licenses=None) -> dict[str, pd.DataFrame]:
        """Same tables, with xref/synonym rows from redistribution-prohibited
        sources/licenses removed (e.g. DrugBank). Compound NODES are kept; only
        the attributions that can't be shared are stripped. A compound that
        exists only because of a prohibited source will keep its node but lose
        the offending rows - check whether any node is left with no shareable
        evidence and decide how to handle that for the public release.
        """
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

    # ---- stats ------------------------------------------------------
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


# ----------------------------------------------------------------------
# Demo of the intended build pattern (replace with real source loaders).
# Run:  python hub.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    hub = CompoundHub()

    # Phase 2: a structure-bearing record (e.g. from Repurposing Hub samples)
    sid = hub.add_structure(
        smiles="CC(=O)Oc1ccccc1C(=O)O",      # aspirin
        source="repurposing_hub_samples",
        preferred_name="aspirin",
    )
    # Phase 3: its external ids - BRD is just one of them, no special status
    hub.add_xref(sid, "BRD", "BRD-K12345678", source="repurposing_hub_samples", license="custom")
    hub.add_xref(sid, "PUBCHEM_CID", "2244", source="repurposing_hub_samples", license="custom")
    # Phase 4 direct: names attached to the structure record
    hub.add_synonym(sid, "Aspirin", source="repurposing_hub_samples", license="custom")
    hub.add_synonym(sid, "acetylsalicylic acid", source="repurposing_hub_samples", license="custom")

    # Phase 4 indirect: a synonym-only source links in via a shared xref...
    target = hub.resolve(xrefs=[("PUBCHEM_CID", "2244")])
    assert target == sid
    hub.add_synonym(target, "ASA", source="katdb", license="open")
    # ...or via a shared name
    assert hub.resolve(names=["acetylsalicylic acid"]) == sid

    # the SAME compound arriving as a salt must dedup to the SAME node
    sid_salt = hub.add_structure(smiles="CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
                                 source="other", preferred_name="aspirin sodium")
    assert sid_salt == sid, "salt should dedup to the parent node"

    # a DrugBank-sourced synonym (will be excluded from the public build)
    hub.add_synonym(sid, "2-acetoxybenzoic acid", source="drugbank", license="cc-by-nc")

    # Phase 5: a name with no structure and no match is KEPT as an orphan
    orphan = hub.resolve(names=["mystery compound x"])
    if orphan is None:
        orphan = hub.add_orphan(source="some_list", preferred_name="mystery compound x")
        hub.add_synonym(orphan, "mystery compound x", source="some_list", license="open")

    print("aspirin node:", sid, "| salt ->", sid_salt, "| orphan:", orphan)
    print("summary:", hub.summary())

    tables = hub.export()
    for name, df in tables.items():
        print(f"\n## {name}\n{df.to_string(index=False)}")

    # Phase 6: redistributable build drops DrugBank-sourced rows
    redist = hub.export_redistributable(prohibited_sources=["drugbank"])
    full_syn = len(tables["synonyms"])
    redist_syn = len(redist["synonyms"])
    print(f"\nsynonym rows  full={full_syn}  redistributable={redist_syn}  "
          f"(dropped {full_syn - redist_syn} DrugBank-sourced)")
