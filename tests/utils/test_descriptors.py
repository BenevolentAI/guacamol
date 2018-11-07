from rdkit import Chem

from guacamol.utils.descriptors import num_atoms_of_type, num_atoms_of_type_fn, num_atoms


def test_num_atoms():
    smiles = 'CCOC(CCC)'
    mol = Chem.MolFromSmiles(smiles)
    assert num_atoms(mol) == 21


def test_num_atoms_does_not_change_mol_instance():
    smiles = 'CCOC(CCC)'
    mol = Chem.MolFromSmiles(smiles)

    assert mol.GetNumAtoms() == 7
    num_atoms(mol)
    assert mol.GetNumAtoms() == 7


def test_count_c_atoms():
    smiles = 'CCOC(CCC)'
    mol = Chem.MolFromSmiles(smiles)
    assert num_atoms_of_type(mol, 'C') == 6


def test_count_h_atoms():
    smiles = 'CCOC(CCC)'
    mol = Chem.MolFromSmiles(smiles)
    assert num_atoms_of_type(mol, 'H') == 14


def test_count_h_atoms_does_not_change_mol_instance():
    smiles = 'CCOC(CCC)'
    mol = Chem.MolFromSmiles(smiles)

    assert mol.GetNumAtoms() == 7
    num_atoms_of_type(mol, 'H')
    assert mol.GetNumAtoms() == 7


def test_closure():
    closure = num_atoms_of_type_fn(element='C')
    mol1 = Chem.MolFromSmiles('CCOCC')
    mol2 = Chem.MolFromSmiles('CNNC')

    assert closure(mol1) == 4
    assert closure(mol2) == 2
