"""
chembl.py
=========
Phase 4+5 extension: Add synonyms from ChEMBL (via InChI Key lookup).

Two-layer design:
  1. add_chembl(hub, cache_path) — reads a pre-fetched TSV cache and
     attaches synonyms/ChEMBL IDs to existing SynDRA nodes. Fast, offline.
  2. fetch_chembl_cache(inchikeys, cache_path) — hits the ChEMBL REST API
     (one request per IK) and writes/updates the cache TSV.

Cache format (TSV, plain text):
  inchikey <TAB> chembl_id <TAB> synonyms (pipe-separated)

API strategy (one call per IK; ChEMBL allows ~4 req/s):
  GET /api/data/molecule?molecule_structures__standard_inchi_key={ik}&format=json
  → molecule_synonyms list (preferred-case names) + pref_name
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

from hub import CompoundHub

_SOURCE = "chembl"
_LICENSE = "cc-by-sa-4.0"
_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_SLEEP = 0.25  # seconds between requests (~4 req/s, safely under ChEMBL limit)
_MAX_SYN_LEN = 120


def _looks_like_name(s: str) -> bool:
    if not s or len(s) > _MAX_SYN_LEN:
        return False
    if s.startswith("InChI="):
        return False
    return True


def _fetch_molecule(ik: str, timeout: int = 15) -> dict | None:
    url = f"{_BASE}/molecule?molecule_structures__standard_inchi_key={ik}&format=json"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            data = json.loads(resp.read())
        mols = data.get("molecules", [])
        return mols[0] if mols else None
    except urllib.error.HTTPError as e:
        if e.code not in (404, 429):
            raise
        return None
    except Exception:
        return None


def fetch_chembl_cache(inchikeys: list[str], cache_path: str,
                       quiet: bool = False) -> None:
    """Fetch ChEMBL synonyms for inchikeys not yet in cache and append results.

    Only IKs absent from the cache are fetched; existing entries are preserved.
    """
    cache_path = Path(cache_path)

    existing_iks: set[str] = set()
    rows: list[dict] = []
    if cache_path.exists():
        with open(cache_path, encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                existing_iks.add(row["inchikey"])
                rows.append(row)

    unique_iks = list(dict.fromkeys(ik for ik in inchikeys if ik))
    missing = [ik for ik in unique_iks if ik not in existing_iks]

    if not quiet:
        print(f"  ChEMBL fetch: {len(unique_iks):,} IKs requested, "
              f"{len(missing):,} not in cache")
    if not missing:
        return

    new_rows: list[dict] = []
    for i, ik in enumerate(missing):
        mol = _fetch_molecule(ik)
        if mol is not None:
            chembl_id = mol.get("molecule_chembl_id", "") or ""
            pref = (mol.get("pref_name") or "").strip()
            raw_syns = mol.get("molecule_synonyms", []) or []
            names: list[str] = []
            if pref and _looks_like_name(pref):
                names.append(pref)
            for s in raw_syns:
                nm = (s.get("molecule_synonym") or "").strip()
                if nm and nm not in names and _looks_like_name(nm):
                    names.append(nm)
            new_rows.append({
                "inchikey": ik,
                "chembl_id": chembl_id,
                "synonyms": "|".join(names),
            })
        time.sleep(_SLEEP)
        if not quiet and (i + 1) % 100 == 0:
            print(f"    … {i + 1}/{len(missing)} fetched")

    if not quiet:
        print(f"  ChEMBL: fetched {len(new_rows)} molecules "
              f"({len(missing) - len(new_rows)} not found in ChEMBL)")

    rows.extend(new_rows)
    with open(cache_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["inchikey", "chembl_id", "synonyms"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    if not quiet:
        print(f"  ChEMBL cache: {len(rows):,} total entries → {cache_path.name}")


def add_chembl(hub: CompoundHub, cache_path: str) -> int:
    """Read pre-fetched ChEMBL cache and attach synonyms/IDs to hub nodes.

    For IKs already in the hub: enriches with synonyms and ChEMBL ID.
    For IKs not yet in the hub: creates an IK-only node (has_structure=False)
    when there is at least a ChEMBL ID or synonym to attach.
    Returns number of compounds enriched or created.
    """
    cache_path = Path(cache_path)
    if not cache_path.exists():
        print(f"  ChEMBL cache not found, skipping: {cache_path.name}")
        return 0

    csv.field_size_limit(10_000_000)
    n_enriched = 0
    n_created = 0

    with open(cache_path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ik        = (row.get("inchikey") or "").strip()
            chembl_id = (row.get("chembl_id") or "").strip()
            syns_raw  = (row.get("synonyms") or "")
            if not ik:
                continue

            syns = [s.strip() for s in syns_raw.split("|") if s.strip()]

            sid = hub._inchikey_to_id.get(ik)
            if not sid:
                # Create IK-only node when we have useful data to attach
                if not syns and not chembl_id:
                    continue
                pref = syns[0] if syns else None
                sid = hub.add_structure(inchikey=ik, preferred_name=pref,
                                        source=_SOURCE)
                if not sid:
                    continue
                n_created += 1

            for s in syns:
                hub.add_synonym(sid, s, source=_SOURCE, license=_LICENSE)
            if chembl_id:
                hub.add_xref(sid, "CHEMBL_ID", chembl_id,
                             source=_SOURCE, license=_LICENSE)
            n_enriched += 1

    print(f"  ChEMBL: enriched {n_enriched:,} compounds "
          f"({n_created} new IK-only nodes created) from cache")
    return n_enriched
