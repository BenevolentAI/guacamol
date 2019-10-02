import subprocess

from guacamol.scoring_function import MoleculewiseScoringFunction
from rdkit import Chem
from similaritymetrics.scoring.base import LigandScorer
from similaritymetrics.utils import io
from similaritymetrics.utils.conformers import generate_conformations_list


class SimilarityMetricScoringFunction(MoleculewiseScoringFunction):
    def __init__(self, method, options, ref_state, num_conformations):
        self.ligandscorer = LigandScorer(method, scoring_kwargs=options)
        self.num_conformations = num_conformations
        self.ref_mol = io.load_from_file_path(ref_state)
    
    def raw_score(self, smiles: str) -> float:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        mol_confs = generate_conformations_list(mol, self.num_conformations)
        try:
            conf_scores = self.ligandscorer(mol_confs, self.ref_mol)
            # Transform to tuple (list of scores, list of mols) à la Lawrence Phillips's legendary `lod2dol`.
            conf_scores, _ = list(zip(*conf_scores))
        except subprocess.TimeoutExpired:
            conf_scores = 0.0  # zero reward
        return max(conf_scores)