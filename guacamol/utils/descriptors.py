from typing import Callable

from rdkit import Chem
from rdkit.Chem import Descriptors, Mol, rdMolDescriptors


def logP(mol: Mol) -> float:
    return Descriptors.MolLogP(mol)


def qed(mol: Mol) -> float:
    return Descriptors.qed(mol)


def tpsa(mol: Mol) -> float:
    return Descriptors.TPSA(mol)


def bertz(mol: Mol) -> float:
    return Descriptors.BertzCT(mol)


def mol_weight(mol: Mol) -> float:
    return Descriptors.MolWt(mol)


def num_H_donors(mol: Mol) -> int:
    return Descriptors.NumHDonors(mol)


def num_H_acceptors(mol: Mol) -> int:
    return Descriptors.NumHAcceptors(mol)


def num_rotatable_bonds(mol: Mol) -> int:
    return Descriptors.NumRotatableBonds(mol)


def num_aromatic_rings(mol: Mol) -> int:
    return rdMolDescriptors.CalcNumAromaticRings(mol)


def num_atoms(mol: Mol) -> int:
    """
    Returns the total number of atoms, H included
    """
    mol = Chem.AddHs(mol)
    return mol.GetNumAtoms()


def num_atoms_of_type(mol: Mol, element: str) -> int:
    """
    Count the number of atoms of a given type.

    Args:
        mol: molecule
        element: element type to look for, such as 'H' or 'C'

    Returns:
        The number of atoms of the given type.
    """
    # if the molecule contains H atoms, they may be implicit, so add them
    if element == 'H':
        mol = Chem.AddHs(mol)

    return sum(1 for a in mol.GetAtoms() if a.GetSymbol() == element)


def num_atoms_of_type_fn(element: str) -> Callable[[Mol], int]:
    """
    Return the function counting the number of atoms of a given element in molecules (closure).
    """
    def counter(mol: Mol) -> int:
        return num_atoms_of_type(mol, element)
    return counter
