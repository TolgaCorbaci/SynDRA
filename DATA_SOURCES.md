# Data sources

SynDRA integrates five external resources. For a reproducible, citable release,
**fill in the exact release version and download date for each file**, and
**confirm that the license permits redistribution** of the derived merged table
before publishing it. (Editors and reviewers routinely check this.)

| Source | File | URL | Version / release | Downloaded | License / terms |
|--------|------|-----|-------------------|------------|-----------------|
| **Therapeutic Targets Database (TTD)** | `P1-04-Drug_synonyms.txt` | https://idrblab.org/ttd/ | _add (e.g. TTD 2024)_ | _add_ | Free for academic use — **verify redistribution terms** |
| **PRISM Repurposing** | `PRISM_drug_synonyms.csv` | https://github.com/broadinstitute/prism_repurposing | _add_ | _add_ | Broad / DepMap terms — **verify** |
| **LINCS 2020 (CMap)** | `compoundinfo_beta.txt` | https://clue.io/releases/data-dashboard | _add (e.g. 2020 beta)_ | _add_ | clue.io / LINCS data-use terms — **verify (commonly CC-BY)** |
| **DrugCentral** | `drugcentral.dump.*.sql` | https://drugcentral.org | _add_ | _add_ | CC BY-SA 4.0 |
| **DrugBank Open Vocabulary** | `drugbank vocabulary.csv` | https://go.drugbank.com/releases/latest#open-data | 2021 | _add_ | CC BY-NC 4.0 — **excluded from redistributable outputs** |
| **PubChem** | `pubchem_synonyms.tsv` *(generated)* | https://pubchem.ncbi.nlm.nih.gov | PUG REST API | _add_ | CC0 (public domain) |
| **ChEMBL** | `chembl_synonyms.tsv` *(generated)* | https://www.ebi.ac.uk/chembl | REST API (latest) | _add_ | CC BY-SA 4.0 |
| **Drug Repurposing Hub** | `repurposing_hub_annotation.txt` | https://repo-hub.broadinstitute.org | _add (e.g. 2025-08-18)_ | _add_ | Non-commercial use only — cite Corsello *et al.* Nat Med 2017 |

## Reproducing the build

Raw source files belong under `synonyms/input/`. The DrugCentral SQL dump and
Repurposing Hub annotation file are excluded from git (see `.gitignore`) due to
size or license constraints. The PubChem synonym cache (`pubchem_synonyms.tsv`) is generated automatically by
`build_all.py` on first run (requires internet). The ChEMBL synonym cache
(`chembl_synonyms.tsv`) is populated by targeted per-IK lookups; run
`build/chembl.py` directly to expand coverage beyond the pre-seeded entries.

The Repurposing Hub annotation file is optional: place it at
`synonyms/input/repurposing_hub_annotation.txt` and run `make build-web` to
add `clinical_phase`, `disease_area`, and `indication` to the web portal. If
absent, the build skips it gracefully.

## Frozen snapshot

For the paper, deposit the exact input snapshot **and** the released outputs
(`syndra_redistributable_compounds.parquet`, `syndra_redistributable_synonyms.parquet`,
`syndra_redistributable_xrefs.parquet`) to Zenodo to obtain a versioned DOI.
Record that DOI in `README.md` and `CITATION.cff`.
