"""
build_all.py
============
Master build script: Phases 1-6.

Run from the project root:
  python build/build_all.py

Outputs written to outputs/ (parquet + CSV).
"""

from __future__ import annotations

import os
import sys

# Ensure build/ modules are importable regardless of cwd
_BUILD_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT_DIR = os.path.dirname(_BUILD_DIR)
sys.path.insert(0, _BUILD_DIR)

from hub import CompoundHub
import compounds as ph2
import xrefs as ph3
import synonyms_build as ph45
import licensing as ph6

# ---------------------------------------------------------------------------
# Input paths (all relative to project root)
# ---------------------------------------------------------------------------
INPUT = os.path.join(_ROOT_DIR, "synonyms", "input")

LINCS_PATH = os.path.join(INPUT, "compoundinfo_beta.txt")
TTD_PATH = os.path.join(INPUT, "P1-04-Drug_synonyms.txt")
PRISM_PATH = os.path.join(INPUT, "PRISM_drug_synonyms.csv")
KATDB_PATH = os.path.join(INPUT, "L1000_BRD_name_translated_drug_list.csv")
REPHUB_PATH = os.path.join(INPUT, "repurposing_samples.txt")   # optional; skip if absent

OUTPUT_DIR = os.path.join(_ROOT_DIR, "outputs")


def main():
    print("=" * 60)
    print("SynDRA Build Pipeline  (structure-anchored, syndra_id)")
    print("=" * 60)

    hub = CompoundHub()

    # ------------------------------------------------------------------
    # Phase 2: Canonical compound nodes from structure-bearing sources
    # ------------------------------------------------------------------
    print("\n--- Phase 2: Canonical compound nodes ---")
    df_lincs = ph2.load_lincs(hub, LINCS_PATH)
    ph2.load_repurposing_hub_samples(hub, REPHUB_PATH)   # no-op if absent

    # ------------------------------------------------------------------
    # Phase 3a: LINCS cross-references (structure-resolved, so safe first)
    # ------------------------------------------------------------------
    print("\n--- Phase 3: Cross-references ---")
    ph3.add_lincs_xrefs(hub, df_lincs)

    # ------------------------------------------------------------------
    # Phase 4+5: Synonyms and orphan handling
    # (PRISM/TTD xrefs follow AFTER because name-resolution needs synonyms)
    # ------------------------------------------------------------------
    print("\n--- Phase 4+5: Synonyms + orphan handling ---")
    ph45.add_lincs_synonyms(hub, df_lincs)
    ph45.add_katdb_synonyms(hub, KATDB_PATH)

    # Parse TTD (needed for both synonyms and xrefs)
    print("  Parsing TTD …")
    df_ttd = ph45.parse_ttd(TTD_PATH)
    n_ttd_orphans = ph45.add_ttd_synonyms(hub, df_ttd)
    n_prism_orphans = ph45.add_prism_synonyms(hub, PRISM_PATH)

    # ------------------------------------------------------------------
    # Phase 3b: PRISM + TTD xrefs (after synonyms so name-resolve works)
    # ------------------------------------------------------------------
    print("\n--- Phase 3b: PRISM + TTD cross-references (post-synonym) ---")
    ph3.add_prism_xrefs(hub, PRISM_PATH)
    ph3.add_ttd_xrefs(hub, df_ttd)

    total_orphans = n_ttd_orphans + n_prism_orphans
    print(f"\n[Phase 5 summary]  orphan nodes created={total_orphans}"
          " (old pipeline would have dropped these)")

    # ------------------------------------------------------------------
    # Phase 6: License tagging + dual builds
    # ------------------------------------------------------------------
    print("\n--- Phase 6: Dual build outputs ---")
    ph6.write_outputs(hub, OUTPUT_DIR)

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    print("\n=== Build complete ===")
    s = hub.summary()
    for k, v in s.items():
        print(f"  {k}: {v}")
    print(f"\nOutputs written to: {OUTPUT_DIR}")
    print("  syndra_full_*.parquet / .csv")
    print("  syndra_redistributable_*.parquet / .csv")

    if hub.ambiguous_names:
        print(f"\n[AUDIT] {len(hub.ambiguous_names)} synonym(s) resolve to >1 node.")
        print("  Review hub.ambiguous_names for potential merge errors.")

    return hub


if __name__ == "__main__":
    main()
