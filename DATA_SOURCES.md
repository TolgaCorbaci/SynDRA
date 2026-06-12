# Data sources

SynDRA integrates four external resources. For a reproducible, citable release,
**fill in the exact release version and download date for each file**, and
**confirm that the license permits redistribution** of the derived merged table
before publishing it. (Editors and reviewers routinely check this.)

| Source | File | URL | Version / release | Downloaded | License / terms |
|--------|------|-----|-------------------|------------|-----------------|
| **KatDB** (Kat Koler) | `L1000_BRD_name_translated_drug_list.csv` | internal / lab resource | _add_ | _add_ | Lab-internal — obtain permission + attribution before redistribution |
| **Therapeutic Targets Database (TTD)** | `P1-04-Drug_synonyms.txt` | https://idrblab.org/ttd/ | _add (e.g. TTD 2024)_ | _add_ | Free for academic use — **verify redistribution terms** |
| **PRISM Repurposing** | `PRISM_drug_synonyms.csv` | https://github.com/broadinstitute/prism_repurposing | _add_ | _add_ | Broad / DepMap terms — **verify** |
| **LINCS 2020 (CMap)** | `compoundinfo_beta.txt` | https://clue.io/releases/data-dashboard | _add (e.g. 2020 beta)_ | _add_ | clue.io / LINCS data-use terms — **verify (commonly CC-BY)** |

## Reproducing the build

The build notebook (`synonyms/scripts/create_synonym_database.ipynb`) currently
references the raw KatDB, TTD, and PRISM files. These raw inputs are **not all
committed** to the repository. To make the build reproducible, either:

1. **Commit the raw inputs** under `synonyms/input/<source>/` (only if licenses
   permit), or
2. **Document exact retrieval** here so a third party can fetch the same versions,
   and have the notebook read from `synonyms/input/<source>/`.

Replace any absolute paths (e.g. `/Users/...`, `/data/work/tolga/...`) in the
notebooks with repository-relative paths before release.

## Frozen snapshot

For the paper, deposit the exact input snapshot **and** the released outputs
(`merged_200K_drug_synonyms.csv`, `merged_synonyms_parent_augmented.csv`,
`brd_parent_consolidation.csv`) to Zenodo to obtain a versioned DOI. Record that
DOI in `README.md` and `CITATION.cff`.
