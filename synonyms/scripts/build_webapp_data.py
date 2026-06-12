#!/usr/bin/env python3
"""
build_webapp_data.py — generate the web app's data file from the full SynDRA map.

The demo `index.html` ships with an embedded sample so it works standalone. To
serve the complete resource, run this to produce `syndra_data.json`, drop it next
to index.html, and switch the data line in index.html to fetch it (see the comment
there). Output schema matches what the app expects:

    { "compounds": { BRD: {name, target, moa, smiles, inchikey, parent, ttd,
                           pubchem, synonyms[], variants[], n_variants} },
      "index":     { synonym_lowercase: BRD },
      "stats":     {synonyms, broad, parents} }

Requires: pandas, rdkit. Run from synonyms/scripts/.
"""
import argparse
import json
from pathlib import Path

import pandas as pd
from syndra_structural_validation import parent_inchikey

MAX_SYN = 80      # cap synonyms stored per compound (keeps file size sane)
MAX_VARIANTS = 60  # cap displayed salt/form variants per compound


def main(compound_info, synonyms_path, out_path):
    comp = pd.read_csv(compound_info, sep="\t", dtype=str).drop_duplicates(subset=["pert_id"])
    syn = pd.read_csv(synonyms_path, dtype=str)
    syn["synonyms"] = syn["synonyms"].str.lower().str.strip()

    comp["parent"] = comp["canonical_smiles"].map(parent_inchikey)
    brd2parent = dict(zip(comp["pert_id"], comp["parent"]))
    parent2brds = (comp.dropna(subset=["parent"])
                       .groupby("parent")["pert_id"].apply(lambda x: sorted(set(x))).to_dict())
    syn_by_brd = syn.dropna(subset=["synonyms"]).groupby("BROAD_drug_ID")["synonyms"].apply(list).to_dict()
    ttd_by_brd = syn.dropna(subset=["TTD_drug_ID"]).groupby("BROAD_drug_ID")["TTD_drug_ID"].first().to_dict()
    pc_by_brd = syn.dropna(subset=["PubChem_CID"]).groupby("BROAD_drug_ID")["PubChem_CID"].first().to_dict()

    def clean(v):
        return v if (pd.notna(v) and str(v) not in ("", '""', '"-"', "-")) else None

    records, index = {}, {}
    for _, row in comp.iterrows():
        brd = row["pert_id"]
        if pd.isna(brd):
            continue
        pk = brd2parent.get(brd)
        variants = [b for b in parent2brds.get(pk, [brd]) if b != brd] if pk else []
        syns = sorted(set(syn_by_brd.get(brd, [])))
        name = row["cmap_name"] if pd.notna(row["cmap_name"]) else brd
        records[brd] = {
            "name": name, "target": clean(row.get("target")), "moa": clean(row.get("moa")),
            "smiles": clean(row.get("canonical_smiles")), "inchikey": clean(row.get("inchi_key")),
            "parent": pk, "ttd": ttd_by_brd.get(brd), "pubchem": pc_by_brd.get(brd),
            "synonyms": syns[:MAX_SYN], "variants": variants[:MAX_VARIANTS], "n_variants": len(variants),
        }
        for s in syns + [name.lower()]:
            if s:
                index[s] = brd

    data = {"compounds": records, "index": index,
            "stats": {"synonyms": int(syn["synonyms"].nunique()),
                      "broad": int(syn["BROAD_drug_ID"].nunique()),
                      "parents": len(parent2brds)}}
    out_path = Path(out_path)
    out_path.write_text(json.dumps(data, separators=(",", ":")))
    mb = out_path.stat().st_size / 1e6
    print(f"Wrote {out_path}  ({mb:.1f} MB)")
    print(f"  compounds: {len(records):,} | searchable synonyms: {len(index):,}")
    if mb > 25:
        print("  Note: >25 MB. Consider hosting via a release asset or Git LFS, and "
              "showing a loading state on first fetch.")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--compound-info", default="synonyms/scripts/data/compoundinfo_beta.tsv")
    p.add_argument("--synonyms", default="synonyms/scripts/data/merged_200K_drug_synonyms.csv")
    p.add_argument("--out", default="syndra_data.json")
    args = p.parse_args()
    main(args.compound_info, args.synonyms, args.out)
