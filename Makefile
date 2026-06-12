# SynDRA — common tasks. Run from the repo root.
SCRIPTS = synonyms/scripts

.PHONY: setup build benchmark validate enhance reproduce app clean

setup:
	pip install -r requirements.txt

# Rebuild the harmonized synonym map from the raw sources in synonyms/input/
# (regenerates synonyms/scripts/data/merged_200K_drug_synonyms.csv)
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

# Full pipeline: build -> benchmark -> validate -> enhance
reproduce: build benchmark validate enhance

# Launch the Shiny web app locally
app:
	cd $(SCRIPTS) && shiny run app.py

clean:
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.pyc' -delete
