"""
drugbank_vocab.py
=================
Phase 4+5 extension: Add synonyms from DrugBank Open Vocabulary CSV.

File: synonyms/input/Drugbank/drugbank vocabulary.csv
Columns used: DrugBank ID, Common name, Synonyms (pipe-separated), Standard InChI Key

Resolution strategy (mirrors other synonym phases):
  1. If InChI Key present → structure-anchor (add_structure with inchikey only)
  2. Else → name lookup, then orphan

Adds DrugBank IDs as DRUGBANK_ID xrefs (supplements any from DrugCentral).

Run order: after all primary synonym phases so name-resolution benefits from
the full synonym index.
"""

from __future__ import annotations

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from hub import CompoundHub

_SOURCE = "drugbank"
_LICENSE = "cc-by-nc-4.0"


def add_drugbank_vocab(hub: CompoundHub, filepath: str) -> int:
    """Parse DrugBank vocabulary CSV and attach synonyms to hub.

    Returns number of orphan nodes created (no InChI Key + no name match).
    """
    n_added = 0

    with open(filepath, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            did   = (row.get("DrugBank ID") or "").strip()
            name  = (row.get("Common name") or "").strip()
            ik    = (row.get("Standard InChI Key") or "").strip()
            syns_raw = (row.get("Synonyms") or "")

            if not name:
                continue

            all_names = [name] + [s.strip() for s in syns_raw.split("|") if s.strip()]

            # Resolve or create compound node.
            # Without an InChIKey we only merge into existing SynDRA nodes by
            # name — we never create new orphans (would add thousands of
            # biologics/proteins that don't belong in a small-molecule hub).
            if ik:
                sid = hub.add_structure(inchikey=ik, preferred_name=name,
                                        source=_SOURCE)
            else:
                sid = hub.resolve(names=all_names)

            if sid is None:
                continue

            # Attach synonyms
            for n in all_names:
                if n:
                    hub.add_synonym(sid, n, source=_SOURCE, license=_LICENSE)

            # Attach DrugBank ID xref
            if did:
                hub.add_xref(sid, "DRUGBANK_ID", did,
                              source=_SOURCE, license=_LICENSE)

            n_added += 1

    print(f"  DrugBank vocab: {n_added:,} entries processed")
    return 0
