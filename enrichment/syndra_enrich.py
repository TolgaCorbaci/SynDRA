"""
syndra_enrich.py
================
Phase 7: Drug-set enrichment using the DrugEnrichr ORA methodology.

Resolver targets syndra_id (not BROAD_drug_ID).

Key steps:
  1. build_resolver(synonyms_df)  ->  {synonym_norm: syndra_id}
  2. parse_gmt(filepath)          ->  {term: [drug_name, ...]}
  3. harmonize_library(...)       ->  {term: {syndra_id, ...}}
  4. run_ora(query_ids, library, universe_ids)  ->  results DataFrame
  5. enrich(query_names, ...)     ->  top enriched terms
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# Allow running from project root or enrichment/
_HERE = Path(__file__).parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_ROOT / "build"))

from normalize import normalize_name


# ---------------------------------------------------------------------------
# Resolver
# ---------------------------------------------------------------------------

def build_resolver(synonyms_df: pd.DataFrame) -> dict[str, str]:
    """Build {synonym_norm -> syndra_id} lookup from the synonyms table.

    When a synonym is ambiguous (maps to multiple nodes), the first
    alphabetically-sorted syndra_id is used (deterministic; mirrors hub.resolve).
    """
    resolver: dict[str, list[str]] = {}
    for _, row in synonyms_df.iterrows():
        norm = str(row.get("synonym_norm", "")).strip()
        sid = str(row.get("syndra_id", "")).strip()
        if norm and sid:
            resolver.setdefault(norm, []).append(sid)

    return {norm: sorted(sids)[0] for norm, sids in resolver.items()}


def resolve_names(names: list[str], resolver: dict[str, str]) -> tuple[set[str], list[str]]:
    """Resolve a list of drug names to syndra_ids.

    Returns (matched_ids, unresolved_names).
    """
    matched: set[str] = set()
    unresolved: list[str] = []
    for name in names:
        norm = normalize_name(name)
        sid = resolver.get(norm)
        if sid:
            matched.add(sid)
        else:
            unresolved.append(name)
    return matched, unresolved


# ---------------------------------------------------------------------------
# GMT parser
# ---------------------------------------------------------------------------

def parse_gmt(filepath: str) -> dict[str, list[str]]:
    """Parse an Enrichr-style GMT file.

    Format: term_name \\t description \\t drug1 \\t drug2 \\t ...
    Returns {term_name: [drug_name, ...]}.
    """
    library: dict[str, list[str]] = {}
    with open(filepath, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0].strip()
            # parts[1] is description (often empty)
            drugs = [p.strip() for p in parts[2:] if p.strip()]
            if term and drugs:
                library[term] = drugs
    return library


def parse_all_libraries(db_dir: str) -> dict[str, dict[str, list[str]]]:
    """Parse every .txt file in db_dir as a GMT library.
    Returns {library_name: {term: [drugs]}}.
    """
    libs: dict[str, dict[str, list[str]]] = {}
    for fname in sorted(os.listdir(db_dir)):
        if not fname.endswith(".txt"):
            continue
        lib_name = fname[:-4]
        libs[lib_name] = parse_gmt(os.path.join(db_dir, fname))
    return libs


# ---------------------------------------------------------------------------
# Harmonization
# ---------------------------------------------------------------------------

def harmonize_library(
    library: dict[str, list[str]],
    resolver: dict[str, str],
) -> dict[str, set[str]]:
    """Map every drug name in each library term to its syndra_id.
    Terms with zero resolved drugs are kept (empty set) for coverage tracking.
    """
    harmonized: dict[str, set[str]] = {}
    for term, drugs in library.items():
        ids: set[str] = set()
        for drug in drugs:
            norm = normalize_name(drug)
            sid = resolver.get(norm)
            if sid:
                ids.add(sid)
        harmonized[term] = ids
    return harmonized


def library_coverage(
    library: dict[str, list[str]],
    harmonized: dict[str, set[str]],
) -> dict:
    """Compute coverage statistics for a single library."""
    total_members = sum(len(v) for v in library.values())
    resolved_members = sum(len(v) for v in harmonized.values())
    unresolved = total_members - resolved_members
    universe = set().union(*harmonized.values()) if harmonized else set()

    return {
        "n_terms": len(library),
        "total_members": total_members,
        "resolved_members": resolved_members,
        "unresolved_members": unresolved,
        "match_rate": resolved_members / total_members if total_members else 0.0,
        "universe_size": len(universe),
    }


# ---------------------------------------------------------------------------
# ORA (one-sided Fisher's exact + Benjamini-Hochberg FDR)
# ---------------------------------------------------------------------------

def run_ora(
    query_ids: set[str],
    harmonized_library: dict[str, set[str]],
    universe_ids: set[str],
    min_overlap: int = 1,
) -> pd.DataFrame:
    """Run over-representation analysis (ORA) for a query set.

    For each term: Fisher's exact test (one-sided, greater).
    Background N = universe_ids (all resolved drug IDs across the library).

    Returns DataFrame sorted by p-value with BH-adjusted q-values.
    """
    N = len(universe_ids)
    K = len(query_ids & universe_ids)  # query drugs in universe

    rows = []
    for term, term_ids in harmonized_library.items():
        M = len(term_ids)
        if M == 0:
            continue
        overlap_ids = query_ids & term_ids
        k = len(overlap_ids)
        if k < min_overlap:
            continue

        # 2x2 contingency table for one-sided Fisher's exact
        # [[k, K-k], [M-k, N-M-(K-k)]]
        a = k
        b = K - k
        c = M - k
        d = N - M - b
        if d < 0:
            d = 0  # universe may be smaller than expected; clamp

        _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "term": term,
            "overlap": k,
            "term_size": M,
            "query_size": K,
            "background": N,
            "pvalue": pval,
            "overlap_ids": ";".join(sorted(overlap_ids)),
        })

    if not rows:
        return pd.DataFrame(columns=["term", "overlap", "term_size", "query_size",
                                     "background", "pvalue", "qvalue", "overlap_ids"])

    df = pd.DataFrame(rows).sort_values("pvalue").reset_index(drop=True)
    df["qvalue"] = _bh_fdr(df["pvalue"].values)
    return df


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction (standard implementation)."""
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    qvals = pvals * n / ranks
    # Enforce monotonicity from right to left
    qvals_sorted = qvals[order]
    for i in range(n - 2, -1, -1):
        qvals_sorted[i] = min(qvals_sorted[i], qvals_sorted[i + 1])
    result = np.empty(n)
    result[order] = qvals_sorted
    return np.minimum(result, 1.0)


# ---------------------------------------------------------------------------
# High-level entry point
# ---------------------------------------------------------------------------

def enrich(
    query_names: list[str],
    synonyms_df: pd.DataFrame,
    db_dir: str,
    min_overlap: int = 1,
    fdr_threshold: float = 0.05,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], list[str]]:
    """Full enrichment pipeline.

    Returns:
      coverage_df  - per-library coverage statistics
      results      - {library_name: ORA results DataFrame}
      unresolved   - query names that did not resolve to a syndra_id
    """
    resolver = build_resolver(synonyms_df)
    query_ids, unresolved = resolve_names(query_names, resolver)

    all_libs = parse_all_libraries(db_dir)

    coverage_rows = []
    results: dict[str, pd.DataFrame] = {}

    for lib_name, library in all_libs.items():
        harmonized = harmonize_library(library, resolver)
        cov = library_coverage(library, harmonized)
        cov["library"] = lib_name
        coverage_rows.append(cov)

        universe = set().union(*harmonized.values()) if harmonized else set()
        ora = run_ora(query_ids, harmonized, universe, min_overlap=min_overlap)
        results[lib_name] = ora

    coverage_df = pd.DataFrame(coverage_rows)[
        ["library", "n_terms", "total_members", "resolved_members",
         "unresolved_members", "match_rate", "universe_size"]
    ]
    return coverage_df, results, unresolved
