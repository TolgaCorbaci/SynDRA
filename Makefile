# SynDRA — common tasks. Run from the repo root.

.PHONY: setup build-new benchmark-new benchmark-drugbank enrich build-web app-new clean

setup:
	pip install -r requirements.txt

# ── Structure-anchored pipeline (syndra_id canonical key) ──────────────────

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

# Full pipeline: build -> recall benchmark -> enrichment -> web data
reproduce-new: build-new benchmark-new enrich build-web

# Launch the Shiny app (Lookup + Enrichment tabs)
app-new:
	shiny run app/app.py

clean:
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.pyc' -delete
