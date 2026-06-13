# Data sources

SynDRA integrates four external resources. For a reproducible, citable release,
**fill in the exact release version and download date for each file**, and
**confirm that the license permits redistribution** of the derived merged table
before publishing it. (Editors and reviewers routinely check this.)

| Source | File | URL | Version / release | Downloaded | License / terms |
|--------|------|-----|-------------------|------------|-----------------|
| **Therapeutic Targets Database (TTD)** | `P1-04-Drug_synonyms.txt` | https://idrblab.org/ttd/ | _add (e.g. TTD 2024)_ | _add_ | Free for academic use — **verify redistribution terms** |
| **PRISM Repurposing** | `PRISM_drug_synonyms.csv` | https://github.com/broadinstitute/prism_repurposing | _add_ | _add_ | Broad / DepMap terms — **verify** |
| **LINCS 2020 (CMap)** | `compoundinfo_beta.txt` | https://clue.io/releases/data-dashboard | _add (e.g. 2020 beta)_ | _add_ | clue.io / LINCS data-use terms — **verify (commonly CC-BY)** |
| **DrugCentral** | `drugcentral.dump.*.sql` | https://drugcentral.org | _add_ | _add_ | CC BY-SA 4.0 |

## Reproducing the build

Raw source files belong under `synonyms/input/<source>/`. The DrugCentral SQL
dump is excluded from git (listed in `.gitignore`) due to file size. All other
inputs should be placed at the paths shown above before running `make build-new`.

## Frozen snapshot

For the paper, deposit the exact input snapshot **and** the released outputs
(`syndra_redistributable_compounds.parquet`, `syndra_redistributable_synonyms.parquet`,
`syndra_redistributable_xrefs.parquet`) to Zenodo to obtain a versioned DOI.
Record that DOI in `README.md` and `CITATION.cff`.
