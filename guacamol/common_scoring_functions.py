from typing import Callable

from rdkit import Chem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

from guacamol.utils.descriptors import mol_weight, logP, num_H_donors, tpsa, num_atoms, AtomCounter
from guacamol.utils.fingerprints import get_fingerprint
from guacamol.score_modifier import ScoreModifier, MinGaussianModifier, MaxGaussianModifier, GaussianModifier
from guacamol.scoring_function import ScoringFunctionBasedOnRdkitMol, MoleculewiseScoringFunction
from guacamol.utils.chemistry import smiles_to_rdkit_mol, parse_molecular_formula


class RdkitScoringFunction(ScoringFunctionBasedOnRdkitMol):
    """
    Scoring function wrapping RDKit descriptors.
    """

    def __init__(self, descriptor: Callable[[Chem.Mol], float], score_modifier: ScoreModifier = None) -> None:
        """
        Args:
            descriptor: molecular descriptors, such as the ones in descriptors.py
            score_modifier: score modifier
        """
        super().__init__(score_modifier=score_modifier)
        self.descriptor = descriptor

    def score_mol(self, mol: Chem.Mol) -> float:
        return self.descriptor(mol)


class TanimotoScoringFunction(ScoringFunctionBasedOnRdkitMol):
    """
    Scoring function that looks at the fingerprint similarity against a target molecule.
    """

    def __init__(self, target, fp_type, score_modifier: ScoreModifier = None) -> None:
        """
        Args:
            target: target molecule
            fp_type: fingerprint type
            score_modifier: score modifier
        """
        super().__init__(score_modifier=score_modifier)

        self.target = target
        self.fp_type = fp_type
        target_mol = smiles_to_rdkit_mol(target)
        if target_mol is None:
            raise RuntimeError(f'The similarity target {target} is not a valid molecule.')

        self.ref_fp = get_fingerprint(target_mol, self.fp_type)

    def score_mol(self, mol: Chem.Mol) -> float:
        fp = get_fingerprint(mol, self.fp_type)
        return TanimotoSimilarity(fp, self.ref_fp)


class CNS_MPO_ScoringFunction(ScoringFunctionBasedOnRdkitMol):
    """
    CNS MPO scoring function
    """

    def __init__(self, max_logP=6.35584, maxMW=360, min_tpsa=40, max_tpsa=90, max_hbd=0) -> None:
        super().__init__()

        self.logP_gauss = MinGaussianModifier(max_logP, 1)
        self.molW_gauss = MinGaussianModifier(maxMW, 60)
        self.tpsa_maxgauss = MaxGaussianModifier(min_tpsa, 20)
        self.tpsa_mingauss = MinGaussianModifier(max_tpsa, 30)
        self.hbd_gauss = MinGaussianModifier(max_hbd, 2.0)

    def score_mol(self, mol: Chem.Mol) -> float:
        mw = mol_weight(mol)
        lp = logP(mol)
        hbd = num_H_donors(mol)
        mol_tpsa = tpsa(mol)

        o1 = self.tpsa_mingauss(mol_tpsa)
        o2 = self.tpsa_maxgauss(mol_tpsa)
        o3 = self.hbd_gauss(hbd)
        o4 = self.logP_gauss(lp)
        o5 = self.molW_gauss(mw)

        return 0.2 * (o1 + o2 + o3 + o4 + o5)


class IsomerScoringFunction(MoleculewiseScoringFunction):
    """
    Scoring function for closeness to a molecular formula.

    The score penalizes deviations from the required number of atoms for each element type, and for the total
    number of atoms.

    F.i., if the target formula is C2H4, the scoring function is the average of three contributions:
    - number of C atoms with a Gaussian modifier with mu=2, sigma=1
    - number of H atoms with a Gaussian modifier with mu=4, sigma=1
    - total number of atoms with a Gaussian modifier with mu=6, sigma=2
    """

    def __init__(self, molecular_formula: str) -> None:
        """
        Args:
            molecular_formula: target molecular formula
        """
        super().__init__()

        element_occurrences = parse_molecular_formula(molecular_formula)

        total_number_atoms = sum(element_tuple[1] for element_tuple in element_occurrences)

        # scoring functions for each element
        self.functions = [RdkitScoringFunction(descriptor=AtomCounter(element),
                                               score_modifier=GaussianModifier(mu=n_atoms, sigma=1.0))
                          for element, n_atoms in element_occurrences]

        # scoring functions for the total number of atoms
        self.functions.append(RdkitScoringFunction(descriptor=num_atoms,
                                                   score_modifier=GaussianModifier(mu=total_number_atoms, sigma=2.0)))

    def raw_score(self, smiles: str) -> float:
        # return the average of all scoring functions
        return sum(f.score(smiles) for f in self.functions) / len(self.functions)


class SMARTSScoringFunction(ScoringFunctionBasedOnRdkitMol):
    """
    Tests for SMARTS which should be or should not be present in the compound.


    """

    def __init__(self, target: str, inverse=False) -> None:
        """

        :param target: The SMARTS string to match.
        :param inverse: Specifies whether the SMARTS is desired (False) or an antipattern, which we don't want to see
                        in the molecules (inverse=False)
        """
        super().__init__()
        self.inverse = inverse
        self.smarts = target
        self.target = Chem.MolFromSmarts(target)

        assert target is not None

    def score_mol(self, mol: Chem.Mol) -> float:

        matches = mol.GetSubstructMatches(self.target)

        if len(matches) > 0:
            if self.inverse:
                return 0.0
            else:
                return 1.0
        else:
            if self.inverse:
                return 1.0
            else:
                return 0.0
