"""
normalize.py
============
THE single source of truth for drug-name normalization in SynDRA.

Every component that produces or consumes a name-match key MUST import from here:
the synonym-table build, the resolver, and the enrichment harmonizer. If two
components normalize differently, joins silently miss and drugs vanish from
results - which is exactly the failure mode SynDRA exists to fix.

NOTE: reconcile this with the normalization already used in SynDRA_pipeline.ipynb
(the rule that built the current map). The README documents: lowercase, split
multi-synonym strings into separate rows, strip whitespace/formatting. If the
notebook strips additional punctuation or uses a specific split delimiter, port
that here so the existing map and the new build agree. Keep exactly ONE function.
"""

from __future__ import annotations

import re
import unicodedata
from typing import Iterable

# Delimiters seen across SynDRA's synonym sources for multi-synonym cells.
# Confirm per source before trusting; add/remove as needed.
_SYNONYM_DELIMITERS = [";", "|", "\t"]
_DELIM_RE = re.compile("|".join(re.escape(d) for d in _SYNONYM_DELIMITERS))

_WS_RE = re.compile(r"\s+")


def normalize_name(name: str) -> str:
    """Normalize one drug name/synonym to its match key.

    lowercase + Unicode NFKC fold + collapse internal whitespace + strip.
    Keep this minimal and deterministic; do NOT add fuzzy logic here - fuzzy
    matching belongs in a separate, clearly-labeled layer, not the join key.
    """
    if name is None:
        return ""
    s = unicodedata.normalize("NFKC", str(name))
    s = s.strip().lower()
    s = _WS_RE.sub(" ", s)
    return s


def split_synonyms(cell: str) -> list[str]:
    """Split a multi-synonym cell into individual raw names (not yet normalized).

    Splits only on the configured delimiters above. Deliberately does NOT split
    on commas or slashes - many chemical names legitimately contain them
    (e.g. "1,1'-hexamethylenebis[...]"). Adjust per source if needed.
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
    # comma-bearing chemical name must NOT be split
    assert split_synonyms("1,1'-hexamethylenebis[5-(4-chlorophenyl)biguanide]") == \
        ["1,1'-hexamethylenebis[5-(4-chlorophenyl)biguanide]"]
    assert split_synonyms("aspirin; ASA | acetylsalicylic acid") == \
        ["aspirin", "ASA", "acetylsalicylic acid"]
    assert normalize_many(["Aspirin", "aspirin", " ASA "]) == ["aspirin", "asa"]
    print("normalize.py: all checks pass")
