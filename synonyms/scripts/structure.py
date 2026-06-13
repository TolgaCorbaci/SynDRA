"""
structure.py
============
Chemical structure standardization for SynDRA's canonical identity layer.

Every canonical compound node is keyed by a STANDARDIZED InChIKey computed here,
so this module is the foundation the whole hub dedupes on. Getting salt/charge
handling wrong here silently splits or merges compounds downstream, so the
pipeline is deliberately explicit and self-tested (see __main__).

Pipeline per SMILES:
  parse -> normalize functional groups -> take parent fragment (desalt/solvent)
        -> neutralize charges -> InChIKey + skeleton + canonical SMILES + InChI

Default identity = full standardized InChIKey (stereo-aware). The 14-char
skeleton (first block of the key) is returned too, for OPTIONAL salt/stereo-
insensitive merging at the hub level - never merge on skeleton silently.

Tautomer canonicalization is intentionally NOT applied here (known hard problem;
can wrongly merge distinct species). Add it later behind an explicit flag.

Requires: rdkit   (pip install rdkit)
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

# RDKit logs ordinary parse/valence messages to stderr; silence them so a build
# over 100k+ structures doesn't drown the logs. Real failures still surface as
# None returns below.
RDLogger.DisableLog("rdApp.*")

# Reusable components (cheaper to instantiate once).
_NORMALIZER = rdMolStandardize.Normalizer()
_UNCHARGER = rdMolStandardize.Uncharger()


@dataclass(frozen=True)
class StdStructure:
    inchikey: str           # full standardized InChIKey (27 char) - the identity key
    inchikey_skeleton: str  # first 14-char block (connectivity layer)
    standard_inchi: str
    canonical_smiles: str    # desalted / neutralized / canonical


def standardize(smiles: str) -> Optional[StdStructure]:
    """Standardize a SMILES string; return None if it can't be parsed/processed.

    Steps: normalize -> parent fragment (remove salts & solvents) -> uncharge.
    """
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
        mol = rdMolStandardize.FragmentParent(mol)   # strip salt/solvent fragments
        if mol is None or mol.GetNumAtoms() == 0:
            return None
        mol = _UNCHARGER.uncharge(mol)               # neutralize where sensible

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


# ----------------------------------------------------------------------
# Self-test: the salt/charge guarantee the rest of the hub relies on.
# Run:  python structure.py
# ----------------------------------------------------------------------
if __name__ == "__main__":
    parent = "CN(C)CCOC(c1ccccc1)c1ccccc1"          # diphenhydramine
    salt = "CN(C)CCOC(c1ccccc1)c1ccccc1.Cl"         # diphenhydramine hydrochloride
    charged = "C[NH+](C)CCOC(c1ccccc1)c1ccccc1"     # protonated form

    sp, ss, sc = standardize(parent), standardize(salt), standardize(charged)

    print("parent  :", sp.inchikey, "|", sp.canonical_smiles)
    print("salt    :", ss.inchikey, "|", ss.canonical_smiles)
    print("charged :", sc.inchikey, "|", sc.canonical_smiles)

    assert sp and ss and sc
    assert ss.inchikey == sp.inchikey, "desalting failed: salt != parent"
    assert sc.inchikey == sp.inchikey, "uncharging failed: charged != parent"
    assert sp.inchikey_skeleton == sp.inchikey.split("-")[0]

    # garbage in -> None, not an exception
    assert standardize("not_a_smiles") is None
    assert standardize("") is None
    assert standardize(None) is None

    print("\nstructure.py: salt, charge, and failure handling all pass")
