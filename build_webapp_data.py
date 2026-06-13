#!/usr/bin/env python3
"""
build_webapp_data.py — generate syndra_data.json from the new structure-anchored build.

Reads from outputs/syndra_redistributable_*.parquet (produced by build/build_all.py).
Optionally enriches with target/MOA from LINCS compoundinfo_beta.txt.

Output schema:
  {
    "compounds": {
      "SYN0000001": {
        "name": str, "inchikey": str|null, "skeleton": str|null, "smiles": str|null,
        "has_structure": bool, "brd": str|null, "pubchem": str|null, "ttd": str|null,
        "target": str|null, "moa": str|null,
        "synonyms": [str],       # up to MAX_SYN
        "variants": [syndra_id], # other IDs sharing the same skeleton, up to MAX_VARIANTS
        "n_variants": int
      }
    },
    "index":  { synonym_norm: syndra_id },
    "stats":  { synonyms, compounds, with_structure, orphans }
  }

Run from the project root:
  python build_webapp_data.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent
OUTPUTS = ROOT / "outputs"
LINCS_PATH = ROOT / "synonyms" / "input" / "compoundinfo_beta.txt"
OUT_PATH = ROOT / "syndra_data.json"

MAX_SYN = 40
MAX_VARIANTS = 60

_SENTINELS = {"", "nan", "none", "null", "n/a", "na"}


def _clean(v):
    s = str(v).strip() if v is not None else ""
    return None if s.lower() in _SENTINELS else s


def _load(stem):
    p = OUTPUTS / f"{stem}.parquet"
    if p.exists():
        return pd.read_parquet(p)
    c = OUTPUTS / f"{stem}.csv"
    if c.exists():
        return pd.read_csv(c, dtype=str).fillna("")
    raise FileNotFoundError(f"Not found: {p} or {c}. Run build/build_all.py first.")


def main():
    print("Loading SynDRA outputs …")
    compounds_df = _load("syndra_redistributable_compounds")
    xrefs_df     = _load("syndra_redistributable_xrefs")
    synonyms_df  = _load("syndra_redistributable_synonyms")

    # ------------------------------------------------------------------
    # Optional: LINCS target/MOA enrichment (keyed by BRD ID)
    # ------------------------------------------------------------------
    brd_to_meta = {}
    if LINCS_PATH.exists():
        lincs = pd.read_csv(LINCS_PATH, sep="\t", dtype=str).fillna("")
        for _, row in lincs.iterrows():
            brd = _clean(row.get("pert_id", ""))
            if brd:
                brd_to_meta[brd] = {
                    "target": _clean(row.get("target", "")),
                    "moa":    _clean(row.get("moa", "")),
                }
        print(f"  LINCS target/MOA loaded: {len(brd_to_meta):,} entries")

    # ------------------------------------------------------------------
    # Build per-compound xref lookup
    # ------------------------------------------------------------------
    xref_by_sid: dict[str, dict[str, list[str]]] = {}
    for _, row in xrefs_df.iterrows():
        sid  = _clean(str(row.get("syndra_id", "")))
        itype = _clean(str(row.get("id_type", "")))
        val  = _clean(str(row.get("id_value", "")))
        if sid and itype and val:
            xref_by_sid.setdefault(sid, {}).setdefault(itype, [])
            if val not in xref_by_sid[sid][itype]:
                xref_by_sid[sid][itype].append(val)

    # ------------------------------------------------------------------
    # Build per-compound synonym list
    # ------------------------------------------------------------------
    syns_by_sid: dict[str, list[str]] = {}
    for _, row in synonyms_df.iterrows():
        sid = _clean(str(row.get("syndra_id", "")))
        raw = _clean(str(row.get("synonym_raw", "")))
        if sid and raw:
            syns_by_sid.setdefault(sid, [])
            if raw not in syns_by_sid[sid]:
                syns_by_sid[sid].append(raw)

    # ------------------------------------------------------------------
    # Group by skeleton → variants
    # ------------------------------------------------------------------
    skeleton_groups: dict[str, list[str]] = {}
    for _, row in compounds_df.iterrows():
        sid = _clean(str(row.get("syndra_id", "")))
        sk  = _clean(str(row.get("inchikey_skeleton", "")))
        if sid and sk:
            skeleton_groups.setdefault(sk, [])
            skeleton_groups[sk].append(sid)

    # ------------------------------------------------------------------
    # Build records
    # ------------------------------------------------------------------
    print("Building compound records …")
    records: dict[str, dict] = {}
    index:   dict[str, str]  = {}

    for _, row in compounds_df.iterrows():
        sid = _clean(str(row.get("syndra_id", "")))
        if not sid:
            continue

        ik    = _clean(str(row.get("inchikey", "")))
        sk    = _clean(str(row.get("inchikey_skeleton", "")))
        smi   = _clean(str(row.get("canonical_smiles", "")))
        name  = _clean(str(row.get("preferred_name", ""))) or sid
        hs    = str(row.get("has_structure", "false")).lower() == "true"

        xrefs = xref_by_sid.get(sid, {})
        brd   = xrefs.get("BRD", [None])[0]
        pubch = xrefs.get("PUBCHEM_CID", [None])[0]
        ttd   = xrefs.get("TTD", [None])[0]

        # Enrich with LINCS target/MOA if we have a BRD xref
        meta = brd_to_meta.get(brd, {}) if brd else {}
        target = meta.get("target")
        moa    = meta.get("moa")

        # Variants: other syndra_ids sharing the same skeleton
        siblings = skeleton_groups.get(sk, []) if sk else []
        variants = [s for s in siblings if s != sid][:MAX_VARIANTS]

        syns = syns_by_sid.get(sid, [])[:MAX_SYN]

        records[sid] = {
            "name": name, "inchikey": ik, "skeleton": sk, "smiles": smi,
            "has_structure": hs, "brd": brd, "pubchem": pubch, "ttd": ttd,
            "target": target, "moa": moa,
            "synonyms": syns, "variants": variants,
            "n_variants": len(siblings) - 1,   # excludes self
        }

        # Build search index from all synonyms
        for _, srow in synonyms_df[synonyms_df["syndra_id"] == sid].iterrows():
            norm = _clean(str(srow.get("synonym_norm", "")))
            if norm and norm not in index:
                index[norm] = sid

        # Also index BRD IDs (so BRD-K... searches work)
        for brd_val in xrefs.get("BRD", []):
            if brd_val:
                bv = brd_val.lower()
                if bv not in index:
                    index[bv] = sid

    # Rebuild index from normalized synonyms for performance (avoids per-sid iterrows)
    print("  Building search index …")
    index = {}
    for _, row in synonyms_df.iterrows():
        sid  = _clean(str(row.get("syndra_id", "")))
        norm = _clean(str(row.get("synonym_norm", "")))
        if sid and norm and norm not in index:
            index[norm] = sid

    # Also add BRD IDs to the index
    for _, row in xrefs_df.iterrows():
        sid   = _clean(str(row.get("syndra_id", "")))
        itype = _clean(str(row.get("id_type", "")))
        val   = _clean(str(row.get("id_value", "")))
        if sid and itype == "BRD" and val:
            bv = val.lower()
            if bv not in index:
                index[bv] = sid

    # Stats
    n_compounds = len(compounds_df)
    n_with_struct = int(compounds_df["has_structure"].astype(str).str.lower().eq("true").sum())
    n_orphans = n_compounds - n_with_struct
    n_synonyms = len(synonyms_df)

    data = {
        "compounds": records,
        "index": index,
        "stats": {
            "synonyms": n_synonyms,
            "compounds": n_compounds,
            "with_structure": n_with_struct,
            "orphans": n_orphans,
        },
    }

    OUT_PATH.write_text(json.dumps(data, separators=(",", ":")))
    mb = OUT_PATH.stat().st_size / 1e6
    print(f"\nWrote {OUT_PATH}  ({mb:.1f} MB)")
    print(f"  compounds: {len(records):,}")
    print(f"  index keys: {len(index):,}")
    print(f"  stats: synonyms={n_synonyms:,}  compounds={n_compounds:,}"
          f"  with_structure={n_with_struct:,}  orphans={n_orphans:,}")
    if mb > 30:
        print("  Note: >30 MB. Consider hosting via a release asset or Git LFS.")


if __name__ == "__main__":
    main()
