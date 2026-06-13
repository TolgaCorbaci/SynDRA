"""
licensing.py
============
Phase 6: Emit dual build outputs.

  syndra_full.*           - all rows, internal use only
  syndra_redistributable.* - DrugBank-sourced rows excluded, CC-BY-SA

Every xrefs/synonyms row carries source + license (set during Phases 2-5).
This module just filters and writes; the tagging is done at row-insertion time.
"""

from __future__ import annotations

import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

from hub import CompoundHub

# Sources/licenses that cannot be redistributed.
# DrugBank academic license forbids redistribution.
_PROHIBITED_SOURCES = ["drugbank"]
_PROHIBITED_LICENSES = ["cc-by-nc", "cc-by-nc-4.0", "drugbank_academic"]


def write_outputs(hub: CompoundHub, output_dir: str) -> None:
    """Write syndra_full and syndra_redistributable to output_dir as parquet + CSV."""
    os.makedirs(output_dir, exist_ok=True)

    tables_full = hub.export()
    tables_redist = hub.export_redistributable(
        prohibited_sources=_PROHIBITED_SOURCES,
        prohibited_licenses=_PROHIBITED_LICENSES,
    )

    for build_name, tables in [("syndra_full", tables_full),
                                ("syndra_redistributable", tables_redist)]:
        for tbl_name, df in tables.items():
            base = os.path.join(output_dir, f"{build_name}_{tbl_name}")
            _write(df, base)

    _print_summary(tables_full, tables_redist)


def _write(df, base_path: str) -> None:
    """Write DataFrame as parquet (preferred) with CSV fallback."""
    try:
        df.to_parquet(base_path + ".parquet", index=False)
    except Exception:
        pass  # parquet optional; CSV always written
    df.to_csv(base_path + ".csv", index=False)


def _print_summary(full: dict, redist: dict) -> None:
    print("\n=== Phase 6: Licensing / Dual Build ===")
    for tbl in ("compounds", "xrefs", "synonyms"):
        n_full = len(full.get(tbl, []))
        n_redist = len(redist.get(tbl, []))
        dropped = n_full - n_redist
        print(f"  {tbl:<12}  full={n_full:>7}  redistributable={n_redist:>7}"
              f"  dropped={dropped}")

    print("\nRelease license note:")
    print("  syndra_redistributable inherits CC-BY-SA from ChEMBL / DrugCentral sources.")
    print("  DrugBank data excluded from all redistributable outputs.")
    print("  Code: MIT.")
