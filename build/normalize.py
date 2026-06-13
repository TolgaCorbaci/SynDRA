"""
normalize.py
============
THE single source of truth for drug-name normalization in SynDRA.

Every component that produces or consumes a name-match key MUST import from here:
the synonym-table build, the resolver, and the enrichment harmonizer. If two
components normalize differently, joins silently miss and drugs vanish from
results - which is exactly the failure mode SynDRA exists to fix.
"""

from __future__ import annotations

import re
import unicodedata
from typing import Iterable

_SYNONYM_DELIMITERS = [";", "|", "\t"]
_DELIM_RE = re.compile("|".join(re.escape(d) for d in _SYNONYM_DELIMITERS))
_WS_RE = re.compile(r"\s+")


def normalize_name(name: str) -> str:
    """Normalize one drug name/synonym to its match key.

    lowercase + Unicode NFKC fold + collapse internal whitespace + strip.
    Deliberately minimal and deterministic; no fuzzy logic.
    """
    if name is None:
        return ""
    s = unicodedata.normalize("NFKC", str(name))
    s = s.strip().lower()
    s = _WS_RE.sub(" ", s)
    return s


def split_synonyms(cell: str) -> list[str]:
    """Split a multi-synonym cell into individual raw names (not yet normalized).

    Splits only on the configured delimiters. Deliberately does NOT split on
    commas or slashes - many chemical names legitimately contain them.
    """
    if cell is None:
        return []
    parts = _DELIM_RE.split(str(cell))
    return [p.strip() for p in parts if p and p.strip()]


def normalize_many(names: Iterable[str]) -> list[str]:
    """Normalize a list, dropping blanks and deduping while preserving order."""
    out, seen = [], set()
    for n in names:
        k = normalize_name(n)
        if k and k not in seen:
            seen.add(k)
            out.append(k)
    return out


if __name__ == "__main__":
    assert normalize_name("  Aspirin ") == "aspirin"
    assert normalize_name("ACETYLSALICYLIC   ACID") == "acetylsalicylic acid"
    assert normalize_name(None) == ""
    assert split_synonyms("1,1'-hexamethylenebis[5-(4-chlorophenyl)biguanide]") == \
        ["1,1'-hexamethylenebis[5-(4-chlorophenyl)biguanide]"]
    assert split_synonyms("aspirin; ASA | acetylsalicylic acid") == \
        ["aspirin", "ASA", "acetylsalicylic acid"]
    assert normalize_many(["Aspirin", "aspirin", " ASA "]) == ["aspirin", "asa"]
    print("normalize.py: all checks pass")
