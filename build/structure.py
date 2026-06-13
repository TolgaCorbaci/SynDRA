"""
structure.py
============
Chemical structure standardization for SynDRA's canonical identity layer.

Pipeline: parse -> normalize -> parent fragment (desalt) -> uncharge
         -> InChIKey + skeleton + canonical SMILES + InChI

Default identity = full standardized InChIKey (stereo-aware).
Skeleton (first 14 chars) returned for optional salt/stereo-insensitive merging only.
Tautomer canonicalization intentionally NOT applied (known hard problem).

Requires: rdkit
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

_NORMALIZER = rdMolStandardize.Normalizer()
_UNCHARGER = rdMolStandardize.Uncharger()


@dataclass(frozen=True)
class StdStructure:
    inchikey: str
    inchikey_skeleton: str  # first 14-char block (connectivity layer)
    standard_inchi: str
    canonical_smiles: str


def standardize(smiles: str) -> Optional[StdStructure]:
    """Standardize a SMILES string; return None if it can't be parsed/processed."""
    if smiles is None:
        return None
    smiles = str(smiles).strip()
    if not smiles:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        mol = _NORMALIZER.normalize(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        if mol is None or mol.GetNumAtoms() == 0:
            return None
        mol = _UNCHARGER.uncharge(mol)

        inchikey = Chem.MolToInchiKey(mol)
        if not inchikey:
            return None
        std_inchi = Chem.MolToInchi(mol)
        canonical_smiles = Chem.MolToSmiles(mol)
    except Exception:
        return None

    return StdStructure(
        inchikey=inchikey,
        inchikey_skeleton=inchikey.split("-")[0],
        standard_inchi=std_inchi or "",
        canonical_smiles=canonical_smiles or "",
    )


def standardize_to_dict(smiles: str) -> Optional[dict]:
    s = standardize(smiles)
    return asdict(s) if s is not None else None


if __name__ == "__main__":
    parent = "CN(C)CCOC(c1ccccc1)c1ccccc1"
    salt = "CN(C)CCOC(c1ccccc1)c1ccccc1.Cl"
    charged = "C[NH+](C)CCOC(c1ccccc1)c1ccccc1"

    sp, ss, sc = standardize(parent), standardize(salt), standardize(charged)
    assert sp and ss and sc
    assert ss.inchikey == sp.inchikey, "desalting failed"
    assert sc.inchikey == sp.inchikey, "uncharging failed"
    assert standardize("not_a_smiles") is None
    assert standardize(None) is None
    print("structure.py: all checks pass")
