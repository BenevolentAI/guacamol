from math import sqrt
from typing import List

import pytest

from guacamol.common_scoring_functions import IsomerScoringFunction, SMARTSScoringFunction
from guacamol.score_modifier import GaussianModifier
from guacamol.scoring_function import BatchScoringFunction, ArithmeticMeanScoringFunction, GeometricMeanScoringFunction


class MockScoringFunction(BatchScoringFunction):
    """
    Mock scoring function that returns values from an array given in the constructor.
    """

    def __init__(self, values: List[float]) -> None:
        super().__init__()
        self.values = values
        self.index = 0

    def raw_score_list(self, smiles_list: List[str]) -> List[float]:
        start = self.index
        self.index += len(smiles_list)
        end = self.index
        return self.values[start:end]


def test_isomer_scoring_function_returns_one_for_correct_molecule():
    c11h24 = IsomerScoringFunction('C11H24')

    # all those smiles fit the formula C11H24
    smiles1 = 'CCCCCCCCCCC'
    smiles2 = 'CC(CCC)CCCCCC'
    smiles3 = 'CCCCC(CC(C)CC)C'

    assert c11h24.score(smiles1) == 1.0
    assert c11h24.score(smiles2) == 1.0
    assert c11h24.score(smiles3) == 1.0


def test_isomer_scoring_function_penalizes_additional_atoms():
    c11h24 = IsomerScoringFunction('C11H24')

    # all those smiles are C11H24O
    smiles1 = 'CCCCCCCCCCCO'
    smiles2 = 'CC(CCC)COCCCCC'
    smiles3 = 'CCCCOC(CC(C)CC)C'

    # the penalty corresponds to a deviation of 1.0 from the gaussian modifier in the total number of atoms
    penalty_tot_num_atoms = 1.0 - GaussianModifier(mu=0, sigma=2)(1.0)
    expected_score = 1.0 - penalty_tot_num_atoms / 3.0

    assert c11h24.score(smiles1) == pytest.approx(expected_score)
    assert c11h24.score(smiles2) == pytest.approx(expected_score)
    assert c11h24.score(smiles3) == pytest.approx(expected_score)


def test_isomer_scoring_function_penalizes_incorrect_number_atoms():
    c11h24 = IsomerScoringFunction('C12H24')

    # all those smiles fit the formula C11H24O
    smiles1 = 'CCCCCCCCOCCC'
    smiles2 = 'CC(CCOC)CCCCCC'
    smiles3 = 'COCCCC(CC(C)CC)C'

    # the penalty corresponds to a deviation of 1.0 from the gaussian modifier in the number of C atoms
    penalty_tot_num_atoms = 1.0 - GaussianModifier(mu=0, sigma=1)(1.0)
    expected_score = 1.0 - penalty_tot_num_atoms / 3.0

    assert c11h24.score(smiles1) == pytest.approx(expected_score)
    assert c11h24.score(smiles2) == pytest.approx(expected_score)
    assert c11h24.score(smiles3) == pytest.approx(expected_score)


def test_smarts_function():
    mol1 = 'COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc2nccc(n2)c3cn(C)c4ccccc34'
    mol2 = 'Cc1c(C)c2OC(C)(COc3ccc(CC4SC(=O)NC4=O)cc3)CCc2c(C)c1O'
    smarts = '[#7;h1]c1ncccn1'

    scofu1 = SMARTSScoringFunction(target=smarts)
    scofu_inv = SMARTSScoringFunction(target=smarts, inverse=True)

    assert scofu1.score(mol1) == 1.0
    assert scofu1.score(mol2) == 0.0
    assert scofu_inv.score(mol1) == 0.0
    assert scofu_inv.score(mol2) == 1.0

    assert scofu1.score_list([mol1])[0] == 1.0
    assert scofu1.score_list([mol2])[0] == 0.0


def test_arithmetic_mean_scoring_function():
    # define a scoring function returning the mean from two mock functions
    # and assert that it returns the correct values.

    weight_1 = 0.4
    weight_2 = 0.6

    mock_values_1 = [0.232, 0.665, 0.0, 1.0, 0.993]
    mock_values_2 = [0.010, 0.335, 0.8, 0.3, 0.847]

    mock_1 = MockScoringFunction(mock_values_1)
    mock_2 = MockScoringFunction(mock_values_2)

    scoring_function = ArithmeticMeanScoringFunction(scoring_functions=[mock_1, mock_2],
                                                     weights=[weight_1, weight_2])

    smiles = ['CC'] * 5

    scores = scoring_function.score_list(smiles)
    expected = [weight_1 * v1 + weight_2 * v2 for v1, v2 in zip(mock_values_1, mock_values_2)]

    assert scores == expected


def test_geometric_mean_scoring_function():
    # define a scoring function returning the geometric mean from two mock functions
    # and assert that it returns the correct values.

    mock_values_1 = [0.232, 0.665, 0.0, 1.0, 0.993]
    mock_values_2 = [0.010, 0.335, 0.8, 0.3, 0.847]

    mock_1 = MockScoringFunction(mock_values_1)
    mock_2 = MockScoringFunction(mock_values_2)

    scoring_function = GeometricMeanScoringFunction(scoring_functions=[mock_1, mock_2])

    smiles = ['CC'] * 5

    scores = scoring_function.score_list(smiles)
    expected = [sqrt(v1 * v2) for v1, v2 in zip(mock_values_1, mock_values_2)]

    assert scores == expected
