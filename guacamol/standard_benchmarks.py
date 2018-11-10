from guacamol.common_scoring_functions import TanimotoScoringFunction, RdkitScoringFunction, CNS_MPO_ScoringFunction, \
    IsomerScoringFunction
from guacamol.distribution_learning_benchmark import DistributionLearningBenchmark, NoveltyBenchmark
from guacamol.frechet_benchmark import FrechetBenchmark
from guacamol.goal_directed_benchmark import GoalDirectedBenchmark
from guacamol.goal_directed_score_contributions import uniform_specification
from guacamol.score_modifier import MinGaussianModifier, MaxGaussianModifier, ClippedScoreModifier, GaussianModifier
from guacamol.scoring_function import MultiPropertyScoringFunction
from guacamol.utils.descriptors import num_rotatable_bonds, num_aromatic_rings, logP, qed, tpsa


def isomers_c11h24() -> GoalDirectedBenchmark:
    """
    Benchmark to try and get all C11H24 molecules there are.
    There should be 159 if one ignores stereochemistry.
    """

    specification = uniform_specification(159)

    return GoalDirectedBenchmark(name='C11H24 benchmark',
                                 objective=IsomerScoringFunction({'C': 11, 'H': 24}),
                                 contribution_specification=specification)


def isomers_c7h8n2o2() -> GoalDirectedBenchmark:
    """
    Benchmark to try and get 100 isomers for C7H8N2O2.
    """

    specification = uniform_specification(100)

    return GoalDirectedBenchmark(name='C7H8N2O2 benchmark',
                                 objective=IsomerScoringFunction({'C': 7, 'H': 8, 'N': 2, 'O': 2}),
                                 contribution_specification=specification)


def isomers_c9h10n2o2pf2cl() -> GoalDirectedBenchmark:
    """
    Benchmark to try and get 100 isomers for C9H10N2O2PF2Cl.
    """

    specification = uniform_specification(100)

    return GoalDirectedBenchmark(name='C9H10N2O2PF2Cl benchmark',
                                 objective=IsomerScoringFunction(
                                     {'C': 9, 'H': 10, 'N': 2, 'O': 2, 'P': 1, 'F': 2, 'Cl': 1}),
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


def tpsa_benchmark(target: float) -> GoalDirectedBenchmark:
    benchmark_name = f'TPSA (target: {target}) benchmark'
    objective = RdkitScoringFunction(descriptor=tpsa,
                                     score_modifier=GaussianModifier(mu=target, sigma=20.0))

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
    return GoalDirectedBenchmark(name='Chemical beauty benchmark',
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
