import subprocess

from guacamol.scoring_function import MoleculewiseScoringFunction
from rdkit import Chem
from similaritymetrics.scoring.base import LigandScorer
from similaritymetrics.utils import io
from similaritymetrics.utils.conformers import generate_conformations_list
from similaritymetrics.utils.hydrogens import generate_protonation_states
from similaritymetrics.scoring.ligand import protonation_states_to_compounds


class SimilarityMetricScoringFunction(MoleculewiseScoringFunction):
    def __init__(self, method, options, ref_state, num_conformations, num_protonation_states):
        self.ligandscorer = LigandScorer(method, scoring_kwargs=options)
        self.num_conformations = num_conformations
        self.ref_mol = io.load_from_file_path(ref_state)
        self.num_protonation_states = num_protonation_states
    
    def raw_score(self, smiles: str) -> float:
        mol = Chem.MolFromSmiles(smiles)
        mol, segments = generate_protonation_states(mol, self.num_protonation_states)

        mol_confs = generate_conformations_list(mol, self.num_conformations)
        mol_confs = protonation_states_to_compounds(mol_confs, segments)
        conf_scores = self.ligandscorer(mol_confs, self.ref_mol)
        # Transform to tuple (list of scores, list of mols) à la Lawrence Phillips's legendary `lod2dol`.
        conf_scores, _ = list(zip(*conf_scores))
        return max(conf_scores)