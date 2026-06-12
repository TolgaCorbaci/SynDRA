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
- Resolves each compound to a standardized **parent InChIKey** (salt/solvent
  stripped, neutralized), so `imatinib mesylate` and `imatinib` map to one drug
- Provides a **structural consolidation map** that links the many BRD IDs
  representing the same parent compound, recovering replicate L1000 signatures that
  exact-name matching would miss

## Key results

Validating the harmonized resource against the canonical structures in LINCS
`compoundinfo_beta` shows the L1000 catalog is highly redundant: **33,515
structurally-resolved BRD IDs collapse to 11,685 unique parent compounds.** About
**75% of catalog entries are salt/form/stereo variants** of a smaller parent set
(up to 32 BRD IDs for a single compound). SynDRA's parent-aware layer lets a single
query recover all of these instead of one.

Synonym→BRD assignment is **fully unambiguous** — every one of the 199,375 synonyms
maps to exactly one BROAD compound (0 collisions).

On a 527-drug benchmark library, SynDRA matches **71.9%** of drugs to a BROAD ID
versus **58.1%** using the LINCS `cmap_name` field alone — a **+13.9 percentage-point**
gain in recall. See [Validation](#validation).

## Installation

```bash
git clone https://github.com/hidelab/SynDRA.git
cd SynDRA
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

RDKit is required for the structural steps (`pip install rdkit`).

## Usage

```bash
make build       # rebuild the synonym map from raw sources in synonyms/input/
make benchmark   # SynDRA vs LINCS-baseline match rate on the 527-drug library
make validate    # offline structure-based validation of the merge
make enhance     # produce parent-augmented table + BRD consolidation map
make reproduce   # all of the above, in order
make app         # run the interactive web app locally
```

Or call scripts directly from `synonyms/scripts/`:

```bash
python build_synonym_db.py
python benchmark_matching.py
python syndra_structural_validation.py
python syndra_enhance.py
```

The reusable structural core is `parent_inchikey(smiles, connectivity_only=True)`
in `syndra_structural_validation.py` — pass `connectivity_only=False` to preserve
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
        ├── build_synonym_db.py             # build pipeline (raw sources -> merged map)
        ├── benchmark_matching.py           # 527-drug benchmark
        ├── syndra_structural_validation.py # structure-based validation
        ├── syndra_enhance.py               # parent-aug + consolidation
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

Run `make validate` and `make benchmark` to regenerate.

| Metric | Value |
|--------|-------|
| Unique synonyms | 199,375 |
| Unique BROAD IDs | 33,821 |
| Unique TTD IDs | 3,127 |
| Unique PubChem CIDs | 1,111 |
| Synonym-level ambiguity (synonym → >1 BRD) | 0 |
| BRD IDs with a resolved parent structure | 33,515 (85% of catalog) |
| Unique parent compounds | 11,685 |
| Catalog entries that are salt/form/stereo duplicates | ~75% |
| Junk synonyms removed (numeric / ≤2 char) | 2,734 |
| Benchmark match rate — LINCS `cmap_name` baseline | 58.1% (306/527) |
| Benchmark match rate — SynDRA | 71.9% (379/527) |
| Benchmark improvement | +13.9 pp |

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
