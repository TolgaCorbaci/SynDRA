# SynDRA — common tasks. Run from the repo root.
SCRIPTS = synonyms/scripts

.PHONY: setup validate enhance app clean

setup:
	pip install -r requirements.txt

# Offline structure-based validation of the synonym->BRD merge
validate:
	cd $(SCRIPTS) && python syndra_structural_validation.py \
		--synonyms data/merged_200K_drug_synonyms.csv \
		--compound-info data/compoundinfo_beta.tsv

# Produce structure-augmented synonym table + BRD consolidation map
enhance:
	cd $(SCRIPTS) && python syndra_enhance.py

# Launch the Shiny web app locally
app:
	cd $(SCRIPTS) && shiny run app.py

clean:
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find . -name '*.pyc' -delete
