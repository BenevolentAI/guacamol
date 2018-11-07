from guacamol.distribution_learning_benchmark import ValidityBenchmark, UniquenessBenchmark, NoveltyBenchmark
from .mock_generator import MockGenerator


def test_validity_does_not_penalize_duplicates():
    generator = MockGenerator(['CCC', 'CCC'])
    benchmark = ValidityBenchmark(number_samples=2)

    assert benchmark.assess_model(generator).score == 1.0


def test_validity_score_is_proportion_of_valid_molecules():
    generator = MockGenerator(['CCC', 'CC(CC)C', 'invalidMolecule'])
    benchmark = ValidityBenchmark(number_samples=3)

    assert benchmark.assess_model(generator).score == 2.0 / 3.0


def test_uniqueness_penalizes_duplicates():
    generator = MockGenerator(['CCC', 'CCC', 'CCC'])
    benchmark = UniquenessBenchmark(number_samples=3)

    assert benchmark.assess_model(generator).score == 1.0 / 3.0


def test_uniqueness_penalizes_duplicates_with_different_smiles_strings():
    generator = MockGenerator(['C(O)C', 'CCO', 'OCC'])
    benchmark = UniquenessBenchmark(number_samples=3)

    assert benchmark.assess_model(generator).score == 1.0 / 3.0


def test_uniqueness_does_not_penalize_invalid_molecules():
    generator = MockGenerator(['C(O)C', 'invalid1', 'invalid2', 'CCC', 'NCCN'])
    benchmark = UniquenessBenchmark(number_samples=3)

    assert benchmark.assess_model(generator).score == 1.0


def test_novelty_score_is_zero_if_no_molecule_is_new():
    molecules = ['CCOCC', 'NNNNONNN', 'C=CC=C']
    generator = MockGenerator(molecules)
    benchmark = NoveltyBenchmark(number_samples=3, training_set=molecules)

    assert benchmark.assess_model(generator).score == 0.0


def test_novelty_score_is_one_if_all_molecules_are_new():
    generator = MockGenerator(['CCOCC', 'NNNNONNN', 'C=CC=C'])
    benchmark = NoveltyBenchmark(number_samples=3, training_set=['CO', 'CC'])

    assert benchmark.assess_model(generator).score == 1.0


def test_novelty_score_does_not_penalize_duplicates():
    generator = MockGenerator(['CCOCC', 'O(CC)CC', 'C=CC=C', 'CC'])
    benchmark = NoveltyBenchmark(number_samples=3, training_set=['CO', 'CC'])

    # Gets 2 out of 3: one of the duplicated molecules is ignored, so the sampled molecules are
    # ['CCOCC', 'C=CC=C', 'CC'],  and 'CC' is not novel
    assert benchmark.assess_model(generator).score == 2.0 / 3.0


def test_novelty_score_penalizes_invalid_molecules():
    generator = MockGenerator(['CCOCC', 'invalid1', 'invalid2', 'CCCC', 'CC'])
    benchmark = NoveltyBenchmark(number_samples=3, training_set=['CO', 'CC'])

    assert benchmark.assess_model(generator).score == 2.0 / 3.0
