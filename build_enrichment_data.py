#!/usr/bin/env python3
"""
build_enrichment_data.py
Generate enrichment_data.json for client-side ORA in the SynDRA web app.

Format:
  {
    "drugs": ["SYN0000001", ...],          # union of all drugs across all libs
    "libraries": {
      "Library_Name": {
        "universe": N,                     # unique drugs in this library
        "terms": {
          "Term Name": [0, 5, 23, ...]     # integer indices into "drugs" array
        }
      }
    }
  }

Run from project root:
  python build_enrichment_data.py
"""

import json
import sys
from pathlib import Path
from unicodedata import normalize as unorm

import pandas as pd

ROOT = Path(__file__).parent
GMT_DIR  = ROOT / "synonyms" / "input" / "enrichment_databases"
SYN_PATH = ROOT / "outputs" / "syndra_redistributable_synonyms.parquet"
OUT_PATH = ROOT / "enrichment_data.json"

MIN_MEMBERS = 2   # drop terms with fewer resolved members

# Large text-mining / raw-signature libraries excluded from the web bundle
# (too many members → file too large; still included in the Python CLI tool)
_WEB_SKIP = {
    "Geneshot_Predicted_Enrichr",
    "Geneshot_Predicted_GeneRIF",
    "Geneshot_Predicted_Tagger",
    "Geneshot_Predicted_from_AutoRIF",
    "Geneshot_Predicted_from_Co-expression",
    "L1000FWD_Signature_Up",
    "L1000FWD_Signature_Down",
}


# ---------------------------------------------------------------------------
# Normalizer (mirrors build/normalize.py — keeps script self-contained)
# ---------------------------------------------------------------------------
def _norm(name: str) -> str:
    return unorm("NFKC", str(name)).strip().lower()


def build_resolver(syn_path: Path) -> dict[str, str]:
    """synonym_norm → syndra_id from redistributable synonyms."""
    df = pd.read_parquet(syn_path) if syn_path.suffix == ".parquet" else pd.read_csv(syn_path, dtype=str).fillna("")
    resolver: dict[str, str] = {}
    for _, row in df.iterrows():
        sid  = str(row.get("syndra_id", "")).strip()
        norm = str(row.get("synonym_norm", "")).strip()
        if sid and norm and norm not in resolver:
            resolver[norm] = sid
    return resolver


def parse_gmt(filepath: str) -> dict[str, list[str]]:
    """Parse a GMT file → {term_name: [member_name, ...]}."""
    lib: dict[str, list[str]] = {}
    with open(filepath, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term, _desc, *members = parts
            term = term.strip()
            members = [m.strip() for m in members if m.strip()]
            if term and members:
                lib[term] = members
    return lib


def harmonize(lib: dict[str, list[str]], resolver: dict[str, str]) -> dict[str, list[str]]:
    """Replace member names with resolved syndra_ids; drop unresolved members."""
    out: dict[str, list[str]] = {}
    for term, members in lib.items():
        ids: list[str] = []
        seen: set[str] = set()
        for m in members:
            sid = resolver.get(_norm(m))
            if sid and sid not in seen:
                ids.append(sid)
                seen.add(sid)
        if len(ids) >= MIN_MEMBERS:
            out[term] = ids
    return out


def main() -> None:
    if not SYN_PATH.exists():
        csv = SYN_PATH.with_suffix(".csv")
        if not csv.exists():
            print(f"ERROR: synonyms not found at {SYN_PATH}. Run build/build_all.py first.")
            sys.exit(1)

    print("Loading resolver …")
    resolver = build_resolver(SYN_PATH)
    print(f"  {len(resolver):,} synonym entries")

    gmt_files = sorted(GMT_DIR.glob("*.txt"))
    if not gmt_files:
        print(f"ERROR: no .txt files in {GMT_DIR}")
        sys.exit(1)
    print(f"Parsing {len(gmt_files)} libraries …")

    # Collect harmonized term→syndra_id_list per library
    all_drug_ids: set[str] = set()
    libraries_raw: dict[str, dict[str, list[str]]] = {}

    for f in gmt_files:
        lib_name = f.stem
        if lib_name in _WEB_SKIP:
            print(f"  {lib_name:<50s}  [skipped — too large for web bundle]")
            continue
        raw = parse_gmt(str(f))
        harm = harmonize(raw, resolver)
        libraries_raw[lib_name] = harm

        # accumulate all member IDs
        for ids in harm.values():
            all_drug_ids.update(ids)

        n_terms_in  = len(raw)
        n_terms_out = len(harm)
        n_members   = sum(len(ids) for ids in harm.values())
        n_unique    = len({sid for ids in harm.values() for sid in ids})
        pct = n_terms_out / n_terms_in * 100 if n_terms_in else 0
        print(f"  {lib_name:<50s}  {n_terms_out:>5}/{n_terms_in:<5} terms  "
              f"{n_members:>7,} members  {n_unique:>5,} unique drugs  ({pct:.0f}%)")

    # Build global drug index (deterministic order)
    drugs: list[str] = sorted(all_drug_ids)
    drug_idx: dict[str, int] = {sid: i for i, sid in enumerate(drugs)}

    # Convert to integer-indexed format
    libraries_out: dict[str, dict] = {}
    for lib_name, harm in libraries_raw.items():
        universe_ids: set[str] = set()
        terms_idx: dict[str, list[int]] = {}
        for term, ids in harm.items():
            idx = [drug_idx[sid] for sid in ids]
            terms_idx[term] = idx
            universe_ids.update(ids)
        libraries_out[lib_name] = {
            "universe": len(universe_ids),
            "terms": terms_idx,
        }

    payload = {"drugs": drugs, "libraries": libraries_out}
    OUT_PATH.write_text(json.dumps(payload, separators=(",", ":")))

    mb  = OUT_PATH.stat().st_size / 1e6
    n_t = sum(len(v["terms"]) for v in libraries_out.values())
    n_m = sum(len(idx) for v in libraries_out.values() for idx in v["terms"].values())
    print(f"\nWrote {OUT_PATH}  ({mb:.1f} MB)")
    print(f"  drugs in any library : {len(drugs):,}")
    print(f"  total terms          : {n_t:,}")
    print(f"  total member entries : {n_m:,}")
    if mb > 20:
        print("  Note: >20 MB — consider setting MIN_MEMBERS=3 to reduce size.")


if __name__ == "__main__":
    main()
