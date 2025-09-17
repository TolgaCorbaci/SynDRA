# SynDRA: Synonym Mapping for Alignment of Repurposing Therapeutics

**SynDRA** is a unified drug synonym mapping system designed to harmonize identifiers across major biomedical resources.  
It bridges gaps between external drug sources and transcriptomic perturbation datasets such as **LINCS/CMap**, improving drug repurposing workflows by resolving inconsistent naming.

## 🌐 Web App

Try SynDRA interactively online: [https://tolgacorbaci.shinyapps.io/syndra/](https://tolgacorbaci.shinyapps.io/syndra/)


![SynDRA Pipeline Overview](https://github.com/hidelab/SynDRA/blob/main/SynDRA%20figure.png)


---

## 📌 Overview

- Integrates synonyms from **KatDB**, **TTD**, **PRISM**, and **LINCS2020**
- Normalizes and deduplicates synonyms into a single mapping
- Links identifiers across **BRD_IDs**, **TTD_IDs**, and **PubChem CIDs**
- Increases match rates for drug repurposing pipelines (e.g., +8.5% for FO5A benchmark set)

**Result:**  
193,113 unique synonyms mapped to 33,858 BRD_IDs, 2,775 TTD_IDs, and 950 PubChem CIDs.

---

## 🔬 Methods

### Data Sources

1. **KatDB Synonyms**  
   - *Source*: Kat Koler  
   - *File*: `L1000_BRD_name_translated_drug_list.csv`  
   - *Purpose*: Expand recognition of BRD compounds in L1000 assays  
   - *Example*: `BRD-K52256627` → “chlorhexidine,” “N-[4-methylpiperazinyl]-…”

2. **Therapeutic Targets Database (TTD)**  
   - *URL*: [TTD Download](https://idrblab.org/ttd/)  
   - *File*: `P1-04-Drug_synonyms.txt`  
   - *Purpose*: Drug–target linked synonyms  
   - *Example*: `D00AAN` → “d00aan,” chemical descriptors

3. **PRISM Drug Synonyms**  
   - *URL*: [PRISM GitHub](https://github.com/broadinstitute/prism_repurposing)  
   - *File*: `PRISM_drug_synonyms.csv`  
   - *Purpose*: Supports MOA enrichment  
   - *Example*: `PubChem_CID 11314340` → “a-674563”

4. **LINCS 2020 Compound Metadata**  
   - *URL*: [Clue.io Data Dashboard](https://clue.io/releases/data-dashboard)  
   - *File*: `compoundinfo_beta.txt`  
   - *Purpose*: Standard compound identifiers for L1000 perturbation studies  

---

### Data Characteristics

- **Initial entry counts**  
  - katdb_df: 13,176  
  - ttd_df: 299,047  
  - prism_df: 112,784  

- **Initial unique identifiers**  
  - BROAD_drug_IDs: 5,539  
  - TTD_drug_IDs: 30,713  
  - PubChem_CIDs: 1,351  

- **LINCS2020 coverage**  
  - Unique BROAD_drug_IDs: 33,613  
  - Synonyms: 34,234  

- **LINCS2020 + KatDB combined**  
  - Unique BROAD_drug_IDs: 33,858  
  - Synonyms: 45,617  

---

### Data Preparation

- Convert all synonyms to lowercase  
- Split multi-synonym strings into separate rows  
- Strip whitespace and formatting artifacts  

---

### Merging Strategy

1. **Synonym explosion** → one synonym per row  
2. **Outer join on synonyms** → maximize matches  
3. **ID propagation** → forward/backward fill with shared synonyms  
4. **Grouping & aggregation** → unique sets of IDs  
5. **Filter** → drop rows without `BROAD_drug_ID`  

---

## 📊 Results

### Final Dataset

- **Rows**: 193,221  
- **Unique synonyms**: 193,113  
- **Identifiers**:  
  - BROAD_drug_IDs: 33,858  
  - TTD_drug_IDs: 2,775  
  - PubChem_CIDs: 950  

- **Missing values (NaNs)**:  
  - BROAD_drug_ID: 0  
  - synonyms: 0  
  - TTD_drug_ID: 45,737  
  - PubChem_CID: 88,237  

- **Duplicate rows**: 0  

---

### Sample Entries

| BROAD_drug_ID | Synonym                     | TTD_drug_ID | PubChem_CID |
|---------------|-----------------------------|-------------|-------------|
| BRD-K52256627 | chlorhexidine               | D0V4GY      | 9552079     |
| BRD-K52256627 | chlorhexidine, combinations | D0V4GY      | 9552079     |
| BRD-K52256627 | 1,1'-hexamethylenebis[...]  | D0V4GY      | 9552079     |

---

### Match Statistics

| Category                        | Initial Entries | Final Entries |
|---------------------------------|-----------------|---------------|
| BROAD_drug_ID                   | 45,617          | 33,858        |
| TTD_drug_ID                     | 299,047         | 2,775         |
| PubChem_CID                     | 112,784         | 950           |
| **Total merged rows**           | —               | 435,530       |
| **Final unique synonym groups** | —               | 193,221       |
| **Matched identifier instances**| —               | 1,145         |

---

## 🚀 Usage

Clone and run the pipeline:

```bash
git clone https://github.com/hidelab/SynDRA.git
cd SynDRA
jupyter notebook SynDRA_pipeline.ipynb
```
---

## 📄 License

This project is licensed under the MIT License.

---

## 📚 Citation

If you use SynDRA in your work, please cite:

> T. Corbaci et. al., *SynDRA: Synonym Mapping for Alignment of Repurposing Therapeutics* (Brazilian Symposium on Bioinformatics BSB, 2025).

