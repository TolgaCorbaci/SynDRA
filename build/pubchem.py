"""
pubchem.py
==========
Phase 4+5 extension: Add synonyms from PubChem (via InChI Key lookup).

Two-layer design:
  1. add_pubchem(hub, cache_path) — reads from a pre-fetched TSV cache and
     attaches synonyms/CIDs to existing SynDRA nodes. Fast, offline.
  2. fetch_pubchem_cache(inchikeys, cache_path) — hits the PubChem REST API
     in batches and writes the cache TSV. Run once per build update.

Cache format (TSV, plain text):
  inchikey <TAB> cid <TAB> synonyms (pipe-separated)

Batch strategy (two calls per batch of 100 IKs, ~3 min for 31K compounds):
  Step 1: inchikeys → CIDs + synonyms
  Step 2: CIDs → canonical InChIKeys  (maps response back to our IKs)
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

_SOURCE = "pubchem"
_LICENSE = "cc0"

_API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_BATCH_SIZE = 100
_SLEEP = 0.25          # seconds between batches (~4 req/s, safely under PubChem limit)
_MAX_SYN_LEN = 120     # drop synonyms longer than this (usually IUPAC/InChI)
_MIN_SYN_LEN = 2


def _looks_like_name(s: str) -> bool:
    """Keep only human-readable drug names; drop InChI strings, CAS numbers."""
    if len(s) < _MIN_SYN_LEN or len(s) > _MAX_SYN_LEN:
        return False
    if s.startswith("InChI="):
        return False
    # CAS numbers: digits-digits-digit (e.g. "138068-37-8") — 3 parts, all numeric
    parts = s.split("-")
    if len(parts) == 3 and all(p.isdigit() for p in parts):
        return False
    return True


def _get(url: str, timeout: int = 30) -> dict | None:
    """GET a PubChem URL; return parsed JSON or None on 404/error."""
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
        return None
    except Exception:
        return None


def fetch_pubchem_cache(inchikeys: list[str], cache_path: str,
                        quiet: bool = False) -> None:
    """Fetch PubChem synonyms for the given InChI Keys and write to cache TSV.

    Uses a two-step batch strategy to avoid one-by-one API calls:
      Batch → {CID: synonyms}  then  Batch → {CID: canonical InChIKey}
    """
    cache_path = Path(cache_path)
    unique_iks = list(dict.fromkeys(ik for ik in inchikeys if ik))

    if not quiet:
        print(f"  PubChem fetch: {len(unique_iks):,} InChI Keys "
              f"in batches of {_BATCH_SIZE} …")

    rows: list[tuple[str, str, str]] = []  # (inchikey, cid, synonyms_pipe)
    errors = 0
    total_batches = (len(unique_iks) + _BATCH_SIZE - 1) // _BATCH_SIZE

    for batch_idx, i in enumerate(range(0, len(unique_iks), _BATCH_SIZE)):
        batch_iks = unique_iks[i:i + _BATCH_SIZE]
        ik_str = ",".join(batch_iks)

        # Step 1: IK batch → CID + synonyms
        syn_data = _get(f"{_API_BASE}/compound/inchikey/{ik_str}/synonyms/JSON")
        time.sleep(_SLEEP)
        if not syn_data:
            continue

        info_list = syn_data.get("InformationList", {}).get("Information", [])
        if not info_list:
            continue

        cid_to_syns: dict[str, list[str]] = {}
        for entry in info_list:
            cid = str(entry.get("CID", ""))
            syns = [s for s in entry.get("Synonym", []) if _looks_like_name(s)]
            if cid and syns:
                cid_to_syns[cid] = syns

        if not cid_to_syns:
            continue

        # Step 2: CID batch → canonical InChIKeys
        cid_str = ",".join(cid_to_syns.keys())
        ik_data = _get(
            f"{_API_BASE}/compound/cid/{cid_str}/property/InChIKey/JSON")
        time.sleep(_SLEEP)

        if ik_data:
            for prop in ik_data.get("PropertyTable", {}).get("Properties", []):
                cid = str(prop.get("CID", ""))
                ik  = prop.get("InChIKey", "")
                if ik and cid in cid_to_syns:
                    rows.append((ik, cid, "|".join(cid_to_syns[cid])))
        else:
            errors += 1

        if not quiet and batch_idx % 50 == 49:
            print(f"    … {batch_idx + 1}/{total_batches} batches done, "
                  f"{len(rows):,} compounds so far")

    with open(cache_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["inchikey", "cid", "synonyms"])
        w.writerows(rows)

    if not quiet:
        print(f"  PubChem cache written: {len(rows):,} compounds, "
              f"{errors} step-2 errors → {cache_path.name}")


def add_pubchem(hub: CompoundHub, cache_path: str) -> int:
    """Read pre-fetched PubChem cache and attach synonyms to existing hub nodes.

    Matches by InChI Key only; never creates new orphan nodes.
    Returns number of compounds enriched.
    """
    cache_path = Path(cache_path)
    if not cache_path.exists():
        print(f"  PubChem cache not found, skipping: {cache_path.name}")
        return 0

    csv.field_size_limit(10_000_000)
    n_enriched = 0
    with open(cache_path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ik       = (row.get("inchikey") or "").strip()
            cid      = (row.get("cid") or "").strip()
            syns_raw = (row.get("synonyms") or "")

            sid = hub._inchikey_to_id.get(ik)
            if not sid:
                continue

            syns = [s.strip() for s in syns_raw.split("|") if s.strip()]
            for s in syns:
                hub.add_synonym(sid, s, source=_SOURCE, license=_LICENSE)
            if cid:
                hub.add_xref(sid, "PUBCHEM_CID", cid,
                             source=_SOURCE, license=_LICENSE)
            n_enriched += 1

    print(f"  PubChem: enriched {n_enriched:,} compounds from cache")
    return n_enriched
