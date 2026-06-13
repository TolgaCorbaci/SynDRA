# SynDRA: Synonym Mapping for Alignment of Repurposing Therapeutics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
<!-- After you mint a Zenodo DOI (see "Citing SynDRA"), replace the line below: -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

**SynDRA** is a structure-aware drug synonym mapping system that harmonizes drug
identifiers across major biomedical resources. It bridges external drug sources
and transcriptomic perturbation datasets (**LINCS/CMap L1000**), increasing match
rates in signature-based drug repurposing by resolving inconsistent naming **and**
collapsing salt, solvate, and stereo variants onto a single parent compound.

## Portal

Try SynDRA interactively: <https://hidelab.github.io/SynDRA/>

![SynDRA pipeline overview](SynDRA%20figure.png)

## Methodology

SynDRA builds a **structure-anchored compound hub** in six phases:

### Phase 2 — Canonical compound nodes

Every compound with a usable SMILES string is standardized using RDKit:
salt fragments are stripped, charges neutralized, and the canonical InChIKey is
computed. Compounds that share the same full InChIKey are collapsed to a single
`syndra_id` node. This means salt forms, solvates, and co-crystals of the same
parent (e.g. `imatinib mesylate` → `imatinib`) receive one identifier rather than
one per formulation.

### Phase 3 — Cross-reference linking

External identifiers — BRD (LINCS/Broad), PubChem CID, ChEMBL, TTD, UNII,
KEGG, DrugCentral, and others — are attached to their resolved `syndra_id` nodes.
Resolution follows a strict priority: full InChIKey match first, then known xref,
then normalized name, so each external ID is anchored to a chemically unique node.

### Phase 4+5 — Synonym integration and orphan handling

Names from all sources are Unicode-normalized (NFKC), lowercased, and
whitespace-stripped before matching. Sources are integrated in dependency order:

| Source | Type | License |
|--------|------|---------|
| LINCS 2020 (CMap) | `cmap_name` + compound aliases | clue.io terms |
| Therapeutic Targets Database (TTD) | drug names + synonyms | academic |
| PRISM Repurposing | drug names + PubChem synonyms | DepMap terms |
| DrugCentral | INN/synonym table + identifier xrefs | CC BY-SA 4.0 |

Compounds that cannot be resolved to a structure (no usable SMILES or InChIKey
across any source) are retained as **orphan nodes** rather than dropped. This
preserves coverage for biologics, mixtures, and early-stage compounds that are not
yet in structure databases.

### Phase 6 — Dual outputs

Two builds are written to `outputs/`:

- **Full build** (`syndra_full_*`) — all integrated sources.
- **Redistributable build** (`syndra_redistributable_*`) — rows from
  sources with permissive open licenses only, inheriting CC BY-SA 4.0 from
  ChEMBL/DrugCentral. DrugBank rows are excluded from this build.

### Enrichment

Compounds are mapped to **21 pharmacological libraries** spanning targets &
mechanisms of action, gene associations, biological pathways, side effects,
transcriptional signatures, ATC classification, and pharmacogenomics. Client-side
enrichment analysis uses **Fisher's exact test** (one-tailed hypergeometric
probability, computed in log-space to avoid overflow) with
**Benjamini-Hochberg FDR** correction.

## Key statistics

| Metric | Value |
|--------|-------|
| Canonical compound nodes | 59,571 |
| Compounds with resolved structure | 31,200 |
| Name-only (orphan) nodes | 28,371 |
| Total synonym entries | 338,730 |
| Compounds covered by enrichment libraries | 8,591 |
| Enrichment libraries | 21 |
| Enrichment terms | 36,996 |
| Benchmark match rate — LINCS `cmap_name` baseline | 58.1% (306/527) |
| Benchmark match rate — SynDRA | 71.9% (379/527) |
| Benchmark improvement | +13.9 pp |

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
make build-new       # build canonical compounds/xrefs/synonyms -> outputs/
make benchmark-new   # recall benchmark: structure-anchored pipeline
make enrich          # enrichment coverage report + statin control
make build-web       # generate web app data files (syndra_data.json + enrichment_data.json)
make reproduce-new   # full pipeline: build -> benchmark -> enrichment -> web data
make app-new         # run the Shiny app locally (Lookup + Enrichment tabs)
```

## Repository structure

```
SynDRA/
├── README.md
├── LICENSE                         # MIT (code)
├── CITATION.cff                    # how to cite
├── DATA_SOURCES.md                 # provenance, versions, and licenses of inputs
├── requirements.txt
├── Makefile
├── index.html                      # static web portal (search + enrichment)
├── build/                          # build pipeline modules
│   ├── build_all.py                # master build script (phases 2-6)
│   ├── hub.py                      # canonical compound hub
│   ├── compounds.py                # phase 2: structure nodes
│   ├── xrefs.py                    # phase 3: cross-references
│   ├── synonyms_build.py           # phase 4+5: synonyms + orphans
│   ├── drugcentral.py              # phase 4+5: DrugCentral integration
│   ├── licensing.py                # phase 6: dual outputs
│   ├── normalize.py                # name normalization
│   └── structure.py                # SMILES/InChIKey standardization (RDKit)
├── build_webapp_data.py            # generate syndra_data.json
├── build_enrichment_data.py        # generate enrichment_data.json
├── benchmark/                      # recall and gold-standard benchmarks
├── enrichment/                     # ORA enrichment pipeline
├── app/                            # Shiny web app
└── synonyms/
    └── input/                      # raw source files (see DATA_SOURCES.md)
```

## Outputs

| File | Description |
|------|-------------|
| `outputs/syndra_redistributable_compounds.parquet` | Canonical compound nodes (syndra_id, InChIKey, SMILES, preferred name) |
| `outputs/syndra_redistributable_synonyms.parquet` | Synonym table (syndra_id → raw name, normalized name, source) |
| `outputs/syndra_redistributable_xrefs.parquet` | Cross-reference table (syndra_id → external ID type + value) |
| `syndra_data.json` | Web app compound index (search + compound cards) |
| `enrichment_data.json` | Client-side enrichment data (21 libraries, integer-indexed) |

## Data sources

All inputs, with download URLs, release versions, dates, and license terms, are
documented in **[DATA_SOURCES.md](DATA_SOURCES.md)**. Please confirm redistribution
terms for your use before reusing the merged table.

## Limitations

- **Stereochemistry.** The default parent key uses the full InChIKey (including
  stereochemistry layers). Enantiomers and diastereomers are kept distinct by
  default, but salts and solvates are collapsed. Adjust `structure.py` if your
  use case requires connectivity-only collapsing.
- **Orphan coverage.** ~48% of nodes lack a resolved structure (biologics,
  mixtures, proprietary compounds). These nodes carry synonyms and xrefs but
  cannot be validated structurally.
- **Enrichment library coverage.** Libraries are pre-built and static. Novel
  compounds added through orphan nodes will not appear in enrichment results
  unless the library GMT files are updated and `build_enrichment_data.py` is rerun.

## Citing SynDRA

If you use SynDRA, please cite the work (see `CITATION.cff`):

> Corbaci, T. *et al.* SynDRA: Synonym Mapping for Alignment of Repurposing
> Therapeutics. *Brazilian Symposium on Bioinformatics (BSB)*, 2025.

## License

Code is released under the [MIT License](LICENSE). Input datasets retain the
licenses of their original providers — see [DATA_SOURCES.md](DATA_SOURCES.md).
