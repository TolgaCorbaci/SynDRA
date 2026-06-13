# SynDRA — common tasks. Run from the repo root.
SCRIPTS = synonyms/scripts

.PHONY: setup build build-new benchmark benchmark-new validate enhance \
        reproduce enrich build-web app app-new clean

setup:
	pip install -r requirements.txt

# ── Legacy pipeline (BRD-centric) ──────────────────────────────────────────

# Rebuild the harmonized synonym map from the raw sources in synonyms/input/
build:
	cd $(SCRIPTS) && python build_synonym_db.py

# Benchmark match rate: SynDRA vs LINCS cmap_name baseline (527-drug library)
benchmark:
	cd $(SCRIPTS) && python benchmark_matching.py

# Offline structure-based validation of the synonym->BROAD merge
validate:
	cd $(SCRIPTS) && python syndra_structural_validation.py \
		--synonyms data/merged_200K_drug_synonyms.csv \
		--compound-info data/compoundinfo_beta.tsv

# Produce structure-augmented synonym table + BRD consolidation map
enhance:
	cd $(SCRIPTS) && python syndra_enhance.py

# Full legacy pipeline
reproduce: build benchmark validate enhance

# Launch the legacy Shiny web app locally
app:
	cd $(SCRIPTS) && shiny run app.py

# ── New structure-anchored pipeline (syndra_id canonical key) ──────────────

# Build canonical compounds/xrefs/synonyms tables -> outputs/
build-new:
	python build/build_all.py

# Recall lift benchmark: naive vs old BRD-centric vs new structure-anchored
benchmark-new:
	python benchmark/recall.py

# DrugBank gold-standard name->InChIKey evaluation (private; not redistributed)
benchmark-drugbank:
	python benchmark/drugbank_gold.py

# Phase 7: enrichment coverage report + statin control
enrich:
	python enrichment/run_all.py

# Generate web app data files (syndra_data.json + enrichment_data.json)
build-web:
	python build_webapp_data.py
	python build_enrichment_data.py

# Full new pipeline: build -> recall benchmark -> enrichment -> web data
reproduce-new: build-new benchmark-new enrich build-web

# Launch the new Shiny app (Lookup + Enrichment tabs)
app-new:
	shiny run app/app.py

clean:
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.pyc' -delete
