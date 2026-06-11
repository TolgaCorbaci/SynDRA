# SynDRA: Synonym Mapping for Alignment of Repurposing Therapeutics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
<!-- After you mint a Zenodo DOI (see "Citing SynDRA"), replace the line below: -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

**SynDRA** is a structure-aware drug synonym mapping system that harmonizes drug
identifiers across major biomedical resources. It bridges external drug sources
and transcriptomic perturbation datasets (**LINCS/CMap L1000**), increasing match
rates in signature-based drug repurposing by resolving inconsistent naming **and**
collapsing salt, solvate, and stereo variants onto a single parent compound.

## Web app

Try SynDRA interactively: <https://tolgacorbaci.shinyapps.io/syndra/>

![SynDRA pipeline overview](SynDRA%20figure.png)

## What SynDRA does

- Integrates drug synonyms from **KatDB**, **TTD**, **PRISM**, and **LINCS 2020**
- Normalizes and deduplicates synonyms into a single mapping across
  **BRD IDs**, **TTD IDs**, and **PubChem CIDs**
- **(new)** Resolves each compound to a standardized **parent InChIKey** (salt/solvent
  stripped, neutralized), so `imatinib mesylate` and `imatinib` map to one drug
- **(new)** Provides a **structural consolidation map** that links the many BRD IDs
  representing the same parent compound, recovering replicate L1000 signatures that
  exact-name matching would miss

## Key result

Validating the harmonized resource against the canonical structures in LINCS
`compoundinfo_beta` shows the L1000 catalog is highly redundant: **33,515
structurally-resolved BRD IDs collapse to 11,685 unique parent compounds.** About
**75% of catalog entries are salt/form/stereo variants** of a smaller parent set
(up to 32 BRD IDs for a single compound). SynDRA's parent-aware layer lets a single
query recover all of these instead of one.

A structure-based audit of the merge shows synonym→BRD assignment is **99.96%
unambiguous** at the synonym level; of the 76 colliding synonyms, ~37 are benign
salt/form variants and ~36 are genuine ambiguities concentrated on compound-code
strings (e.g. `CHIR-99021`, `AZD-2014`). See [Validation](#validation).

## Installation

```bash
git clone https://github.com/hidelab/SynDRA.git
cd SynDRA
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

RDKit is required for the structural steps (`pip install rdkit`).

## Usage

The harmonized resource is provided directly; the scripts reproduce and validate it.

```bash
# 1. Validate the merge against chemical structure (offline; no internet needed)
make validate

# 2. Produce the structure-augmented synonym table + the BRD consolidation map
make enhance

# 3. Run the interactive web app locally
make app
```

Or call the scripts directly from `synonyms/scripts/`:

```bash
python syndra_structural_validation.py --synonyms data/merged_200K_drug_synonyms.csv \
                                       --compound-info data/compoundinfo_beta.tsv
python syndra_enhance.py
```

The reusable core is `parent_inchikey(smiles, connectivity_only=True)` in
`syndra_structural_validation.py` — pass `connectivity_only=False` to preserve
stereochemistry (see [Limitations](#limitations)).

## Repository structure

```
SynDRA/
├── README.md
├── LICENSE                         # MIT (code)
├── CITATION.cff                    # how to cite
├── DATA_SOURCES.md                 # provenance, versions, and licenses of inputs
├── requirements.txt
├── Makefile
├── SynDRA figure.png
└── synonyms/
    ├── input/                      # raw source files (see DATA_SOURCES.md)
    └── scripts/
        ├── create_synonym_database.ipynb   # build pipeline
        ├── synonym_matching.ipynb          # benchmark matching
        ├── syndra_structural_validation.py # structure-based validation (new)
        ├── syndra_enhance.py               # parent-aug + consolidation (new)
        ├── app.py                          # Shiny web app
        └── data/                           # released outputs
```

## Outputs

| File | Description |
|------|-------------|
| `merged_200K_drug_synonyms.csv` | Harmonized synonym → BRD/TTD/PubChem mapping |
| `merged_synonyms_parent_augmented.csv` | Above + `parent_inchikey`, junk synonyms removed, ambiguous synonyms flagged |
| `brd_parent_consolidation.csv` | `parent_inchikey` → set of BRD IDs that are the same parent compound |

## Data sources

All inputs, with download URLs, release versions, dates, and license terms, are
documented in **[DATA_SOURCES.md](DATA_SOURCES.md)**. Please confirm redistribution
terms for your use before reusing the merged table.

## Validation

Run `make validate` to regenerate. On the released resource:

| Metric | Value |
|--------|-------|
| Unique synonyms | 192,109 |
| Unique BRD IDs | 33,858 |
| BRD IDs with a resolved parent structure | 33,515 (85% of catalog) |
| Unique parent compounds | 11,685 |
| Catalog entries that are salt/form/stereo duplicates | ~75% |
| Synonym-level ambiguity (synonym → >1 BRD) | 76 (0.04%) |
| Junk synonyms removed (numeric / ≤2 char) | 2,413 |

## Limitations

- **BROAD-centric.** Rows without a BRD ID are dropped, so standalone TTD/PubChem
  coverage is intentionally reduced in favor of L1000 linkage.
- **Stereochemistry.** The default parent key uses the InChIKey connectivity block,
  which collapses enantiomers and other stereoisomers. This is usually desirable for
  repurposing but not always (e.g. thalidomide enantiomers differ in activity). Use
  `connectivity_only=False` for a stereo-preserving key.
- **Structure coverage.** ~15% of catalog compounds lack a usable SMILES and are not
  structurally validated.

## Citing SynDRA

If you use SynDRA, please cite the work (see `CITATION.cff`):

> Corbaci, T. *et al.* SynDRA: Synonym Mapping for Alignment of Repurposing
> Therapeutics. *Brazilian Symposium on Bioinformatics (BSB)*, 2025.

## License

Code is released under the [MIT License](LICENSE). Input datasets retain the
licenses of their original providers — see [DATA_SOURCES.md](DATA_SOURCES.md).
