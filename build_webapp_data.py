#!/usr/bin/env python3
"""
build_webapp_data.py — generate syndra_data.json from the new structure-anchored build.

Reads from outputs/syndra_redistributable_*.parquet (produced by build/build_all.py).
Optionally enriches with:
  - target/MOA          from LINCS compoundinfo_beta.txt    (keyed by BRD ID)
  - clinical_phase,
    disease_area,
    indication          from Repurposing Hub annotation file (keyed by pert_iname)

Output schema:
  {
    "compounds": {
      "SYN0000001": {
        "name": str, "inchikey": str|null, "skeleton": str|null, "smiles": str|null,
        "has_structure": bool, "brd": str|null, "pubchem": str|null, "ttd": str|null,
        "approval": str|null,
        "clinical_phase": str|null, "disease_area": str|null, "indication": str|null,
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
from unicodedata import normalize as _unorm

import pandas as pd

ROOT = Path(__file__).parent
OUTPUTS = ROOT / "outputs"
LINCS_PATH = ROOT / "synonyms" / "input" / "compoundinfo_beta.txt"
REPHUB_PATH = ROOT / "synonyms" / "input" / "repurposing_hub_annotation.txt"
OUT_PATH = ROOT / "syndra_data.json"

MAX_SYN = 40
MAX_VARIANTS = 60

_SENTINELS = {"", "nan", "none", "null", "n/a", "na"}
_STATUS_LABEL = {"ONP": "Approved · Rx", "OFP": "Approved · Rx/OTC", "OFM": "Approved · OTC"}


def _clean(v):
    s = str(v).strip() if v is not None else ""
    return None if s.lower() in _SENTINELS else s


def _norm(s: str) -> str:
    return _unorm("NFKC", str(s)).strip().lower()


def _load(stem):
    p = OUTPUTS / f"{stem}.parquet"
    if p.exists():
        return pd.read_parquet(p)
    c = OUTPUTS / f"{stem}.csv"
    if c.exists():
        return pd.read_csv(c, dtype=str).fillna("")
    raise FileNotFoundError(f"Not found: {p} or {c}. Run build/build_all.py first.")


def _dedup_pipe(value: str | None) -> str | None:
    """Deduplicate pipe-separated values while preserving order."""
    if not value:
        return value
    seen: set[str] = set()
    parts: list[str] = []
    for p in value.split("|"):
        p = p.strip()
        if p and p not in seen:
            seen.add(p)
            parts.append(p)
    return " | ".join(parts) if parts else None


def _load_rephub(filepath: Path, syn_index: dict[str, str]) -> dict[str, dict]:
    """Parse Repurposing Hub drug annotation file (! header lines, then TSV).

    Columns used: pert_iname, clinical_phase, moa, target, disease_area, indication.
    Returns {syndra_id: {clinical_phase, disease_area, indication, rephub_moa, rephub_target}}.
    """
    sid_to_rephub: dict[str, dict] = {}
    unresolved = 0

    with open(filepath, encoding="utf-8", errors="replace") as fh:
        headers = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!") or not line.strip():
                continue
            parts = line.split("\t")
            if headers is None:
                headers = [h.strip() for h in parts]
                continue
            row = dict(zip(headers, parts))
            pert = _clean(row.get("pert_iname", ""))
            if not pert:
                continue
            sid = syn_index.get(_norm(pert))
            if not sid:
                unresolved += 1
                continue
            if sid not in sid_to_rephub:
                sid_to_rephub[sid] = {
                    "clinical_phase": _clean(row.get("clinical_phase", "")),
                    "disease_area":   _dedup_pipe(_clean(row.get("disease_area", ""))),
                    "indication":     _dedup_pipe(_clean(row.get("indication", ""))),
                    "rephub_moa":     _dedup_pipe(_clean(row.get("moa", ""))),
                    "rephub_target":  _dedup_pipe(_clean(row.get("target", ""))),
                }

    return sid_to_rephub, unresolved


def main():
    print("Loading SynDRA outputs …")
    compounds_df = _load("syndra_redistributable_compounds")
    xrefs_df     = _load("syndra_redistributable_xrefs")
    synonyms_df  = _load("syndra_redistributable_synonyms")

    # ------------------------------------------------------------------
    # Build synonym index early — needed to resolve RepHub pert_inames.
    # PubChem is excluded from the web search index (it adds ~700K catalog
    # numbers that inflate the JSON without improving drug-name search);
    # PubChem synonyms still live in the redistributable parquets.
    # ------------------------------------------------------------------
    _WEB_INDEX_SOURCES = {"lincs_2020", "ttd", "prism", "drugcentral",
                          "drugbank", "repurposing_hub"}
    print("  Building synonym index …")
    syn_index: dict[str, str] = {}
    for _, row in synonyms_df.iterrows():
        src  = str(row.get("source", ""))
        if src not in _WEB_INDEX_SOURCES:
            continue
        sid  = _clean(str(row.get("syndra_id", "")))
        norm = _clean(str(row.get("synonym_norm", "")))
        if sid and norm and norm not in syn_index:
            syn_index[norm] = sid

    # ------------------------------------------------------------------
    # Optional: LINCS target/MOA enrichment (keyed by BRD ID)
    # ------------------------------------------------------------------
    brd_to_meta: dict[str, dict] = {}
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
    # Optional: Repurposing Hub clinical annotation (keyed by pert_iname)
    # ------------------------------------------------------------------
    sid_to_rephub: dict[str, dict] = {}
    if REPHUB_PATH.exists():
        sid_to_rephub, unresolved = _load_rephub(REPHUB_PATH, syn_index)
        print(f"  Repurposing Hub annotation loaded: {len(sid_to_rephub):,} matched "
              f"({unresolved:,} pert_inames unresolved)")
    else:
        print(f"  Repurposing Hub annotation not found, skipping: {REPHUB_PATH.name}")

    # ------------------------------------------------------------------
    # Build per-compound xref lookup
    # ------------------------------------------------------------------
    xref_by_sid: dict[str, dict[str, list[str]]] = {}
    for _, row in xrefs_df.iterrows():
        sid   = _clean(str(row.get("syndra_id", "")))
        itype = _clean(str(row.get("id_type", "")))
        val   = _clean(str(row.get("id_value", "")))
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

    for _, row in compounds_df.iterrows():
        sid = _clean(str(row.get("syndra_id", "")))
        if not sid:
            continue

        ik   = _clean(str(row.get("inchikey", "")))
        sk   = _clean(str(row.get("inchikey_skeleton", "")))
        smi  = _clean(str(row.get("canonical_smiles", "")))
        name = _clean(str(row.get("preferred_name", ""))) or sid
        hs   = str(row.get("has_structure", "false")).lower() == "true"

        xrefs = xref_by_sid.get(sid, {})
        brd   = xrefs.get("BRD", [None])[0]
        pubch = xrefs.get("PUBCHEM_CID", [None])[0]
        ttd   = xrefs.get("TTD", [None])[0]

        # Approval status from DrugCentral
        dc_status = xrefs.get("DC_STATUS", [None])[0]
        approval  = _STATUS_LABEL.get(dc_status) if dc_status else None

        # Repurposing Hub clinical annotation
        rh = sid_to_rephub.get(sid, {})

        # target/MOA: LINCS BRD match has priority; fall back to RepHub
        lincs_meta = brd_to_meta.get(brd, {}) if brd else {}
        target = lincs_meta.get("target") or rh.get("rephub_target")
        moa    = lincs_meta.get("moa")    or rh.get("rephub_moa")

        siblings = skeleton_groups.get(sk, []) if sk else []
        variants = [s for s in siblings if s != sid][:MAX_VARIANTS]
        syns     = syns_by_sid.get(sid, [])[:MAX_SYN]

        rec = {
            "name": name, "has_structure": hs,
            "inchikey": ik, "skeleton": sk, "smiles": smi,
            "brd": brd, "pubchem": pubch, "ttd": ttd,
            "approval":       approval,
            "clinical_phase": rh.get("clinical_phase"),
            "disease_area":   rh.get("disease_area"),
            "indication":     rh.get("indication"),
            "target": target, "moa": moa,
            "synonyms": syns, "variants": variants,
            "n_variants": len(siblings) - 1,
        }
        records[sid] = {k: v for k, v in rec.items() if v is not None}

    # ------------------------------------------------------------------
    # Build output index: synonym norms + BRD IDs
    # ------------------------------------------------------------------
    index = dict(syn_index)  # already built above
    for _, row in xrefs_df.iterrows():
        sid   = _clean(str(row.get("syndra_id", "")))
        itype = _clean(str(row.get("id_type", "")))
        val   = _clean(str(row.get("id_value", "")))
        if sid and itype == "BRD" and val:
            bv = val.lower()
            if bv not in index:
                index[bv] = sid

    # Stats
    n_compounds   = len(compounds_df)
    n_with_struct = int(compounds_df["has_structure"].astype(str).str.lower().eq("true").sum())
    n_orphans     = n_compounds - n_with_struct
    n_synonyms    = len(synonyms_df)

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
