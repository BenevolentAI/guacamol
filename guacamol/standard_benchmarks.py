from guacamol.common_scoring_functions import TanimotoScoringFunction, RdkitScoringFunction, CNS_MPO_ScoringFunction
from guacamol.distribution_learning_benchmark import DistributionLearningBenchmark, NoveltyBenchmark
from guacamol.frechet_benchmark import FrechetBenchmark
from guacamol.goal_directed_benchmark import GoalDirectedBenchmark
from guacamol.goal_directed_score_contributions import uniform_specification
from guacamol.score_modifier import MinGaussianModifier, MaxGaussianModifier, ClippedScoreModifier, GaussianModifier
from guacamol.scoring_function import MultiPropertyScoringFunction
from guacamol.utils.descriptors import num_rotatable_bonds, num_aromatic_rings, logP, qed, num_atoms_of_type_fn, \
    num_atoms


def isomers_c11h24() -> GoalDirectedBenchmark:
    """
    Benchmark to try and get all C11H24 molecules there are.
    There should be 159 if one ignores stereochemistry.
    """

    number_c = RdkitScoringFunction(descriptor=num_atoms_of_type_fn('C'),
                                    score_modifier=GaussianModifier(mu=11, sigma=1))
    number_h = RdkitScoringFunction(descriptor=num_atoms_of_type_fn('H'),
                                    score_modifier=GaussianModifier(mu=24, sigma=1))

    total_number_atoms = RdkitScoringFunction(descriptor=num_atoms,
                                              score_modifier=GaussianModifier(mu=35, sigma=2))

    specification = uniform_specification(159)

    return GoalDirectedBenchmark(name='C11H24 benchmark',
                                 objective=MultiPropertyScoringFunction([number_c, number_h, total_number_atoms]),
                                 contribution_specification=specification)


def isomers_c7h8n2o2() -> GoalDirectedBenchmark:
    """
    Benchmark to try and get 50 isomers for C7H8N2O2.
    """

    number_c = RdkitScoringFunction(descriptor=num_atoms_of_type_fn('C'),
                                    score_modifier=GaussianModifier(mu=7, sigma=1))
    number_h = RdkitScoringFunction(descriptor=num_atoms_of_type_fn('H'),
                                    score_modifier=GaussianModifier(mu=8, sigma=1))
    number_n = RdkitScoringFunction(descriptor=num_atoms_of_type_fn('N'),
                                    score_modifier=GaussianModifier(mu=2, sigma=1))
    number_o = RdkitScoringFunction(descriptor=num_atoms_of_type_fn('O'),
                                    score_modifier=GaussianModifier(mu=2, sigma=1))

    total_number_atoms = RdkitScoringFunction(descriptor=num_atoms,
                                              score_modifier=GaussianModifier(mu=19, sigma=2))

    specification = uniform_specification(50)

    return GoalDirectedBenchmark(name='C7H8N2O2 benchmark',
                                 objective=MultiPropertyScoringFunction([number_c, number_h, number_n, number_o, total_number_atoms]),
                                 contribution_specification=specification)


def cobimetinib() -> GoalDirectedBenchmark:
    smiles = 'OC1(CN(C1)C(=O)C1=C(NC2=C(F)C=C(I)C=C2)C(F)=C(F)C=C1)C1CCCCN1'

    modifier = ClippedScoreModifier(upper_x=0.7)
    os_tf = TanimotoScoringFunction(smiles, fp_type='FCFP4', score_modifier=modifier)
    os_ap = TanimotoScoringFunction(smiles, fp_type='ECFP6',
                                    score_modifier=MinGaussianModifier(mu=0.75, sigma=0.1))

    rot_b = RdkitScoringFunction(descriptor=num_rotatable_bonds,
                                 score_modifier=MinGaussianModifier(mu=3, sigma=1))

    rings = RdkitScoringFunction(descriptor=num_aromatic_rings,
                                 score_modifier=MaxGaussianModifier(mu=3, sigma=1))

    t_cns = MultiPropertyScoringFunction([os_tf, os_ap, rot_b, rings, CNS_MPO_ScoringFunction()])

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name='Cobimetinib benchmark',
                                 objective=t_cns,
                                 contribution_specification=specification)


def similarity_cns_mpo(smiles, molecule_name) -> GoalDirectedBenchmark:
    benchmark_name = f'{molecule_name} benchmark'
    os_tf = TanimotoScoringFunction(smiles, fp_type='FCFP4')
    os_ap = TanimotoScoringFunction(smiles, fp_type='AP')
    anti_fp = TanimotoScoringFunction(smiles, fp_type='ECFP6',
                                      score_modifier=MinGaussianModifier(mu=0.70, sigma=0.1))

    t_cns = MultiPropertyScoringFunction([os_tf, os_ap, anti_fp, CNS_MPO_ScoringFunction()])

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name=benchmark_name,
                                 objective=t_cns,
                                 contribution_specification=specification)


def similarity(smiles: str, name: str) -> GoalDirectedBenchmark:
    threshold = 0.7
    benchmark_name = f'{name} similarity benchmark'

    modifier = ClippedScoreModifier(upper_x=threshold)
    scoring_function = TanimotoScoringFunction(target=smiles, fp_type='ECFP4', score_modifier=modifier)

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name=benchmark_name,
                                 objective=scoring_function,
                                 contribution_specification=specification)


def logP_benchmark(target: float) -> GoalDirectedBenchmark:
    benchmark_name = f'logP (target: {target}) benchmark'
    objective = RdkitScoringFunction(descriptor=logP,
                                     score_modifier=GaussianModifier(mu=target, sigma=1))

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name=benchmark_name,
                                 objective=objective,
                                 contribution_specification=specification)


def cns_mpo() -> GoalDirectedBenchmark:
    specification = uniform_specification(1, 10, 100)
    return GoalDirectedBenchmark(name='CNS MPO benchmark', objective=CNS_MPO_ScoringFunction(),
                                 contribution_specification=specification)


def qed_benchmark() -> GoalDirectedBenchmark:
    specification = uniform_specification(1, 10, 100)
    return GoalDirectedBenchmark(name='Synthetic accessibility benchmark',
                                 objective=RdkitScoringFunction(descriptor=qed),
                                 contribution_specification=specification)


def median_molecule() -> GoalDirectedBenchmark:
    t_camphor = TanimotoScoringFunction('CC1(C)C2CCC1(C)C(=O)C2', fp_type='ECFP4')
    t_menthol = TanimotoScoringFunction('CC(C)C1CCC(C)CC1O', fp_type='ECFP4')
    median = MultiPropertyScoringFunction([t_menthol, t_camphor])

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name='Median molecules benchmark',
                                 objective=median,
                                 contribution_specification=specification)


def novelty_benchmark(training_set_file: str, number_samples: int) -> DistributionLearningBenchmark:
    smiles_list = [s.strip() for s in open(training_set_file).readlines()]
    return NoveltyBenchmark(number_samples=number_samples, training_set=smiles_list)


def frechet_benchmark(training_set_file: str) -> DistributionLearningBenchmark:
    smiles_list = [s.strip() for s in open(training_set_file).readlines()]
    return FrechetBenchmark(training_set=smiles_list)
