from guacamol.common_scoring_functions import TanimotoScoringFunction, RdkitScoringFunction, CNS_MPO_ScoringFunction, \
    IsomerScoringFunction
from guacamol.distribution_learning_benchmark import DistributionLearningBenchmark, NoveltyBenchmark, KLDivBenchmark
from guacamol.frechet_benchmark import FrechetBenchmark
from guacamol.goal_directed_benchmark import GoalDirectedBenchmark
from guacamol.goal_directed_score_contributions import uniform_specification
from guacamol.score_modifier import MinGaussianModifier, MaxGaussianModifier, ClippedScoreModifier, GaussianModifier
from guacamol.scoring_function import ArithmeticMeanScoringFunction
from guacamol.utils.descriptors import num_rotatable_bonds, num_aromatic_rings, logP, qed, tpsa, bertz, mol_weight, \
    AtomCounter


def isomers_c11h24() -> GoalDirectedBenchmark:
    """
    Benchmark to try and get all C11H24 molecules there are.
    There should be 159 if one ignores stereochemistry.
    """

    specification = uniform_specification(159)

    return GoalDirectedBenchmark(name='C11H24',
                                 objective=IsomerScoringFunction('C11H24'),
                                 contribution_specification=specification)


def isomers_c7h8n2o2() -> GoalDirectedBenchmark:
    """
    Benchmark to try and get 100 isomers for C7H8N2O2.
    """

    specification = uniform_specification(100)

    return GoalDirectedBenchmark(name='C7H8N2O2',
                                 objective=IsomerScoringFunction('C7H8N2O2'),
                                 contribution_specification=specification)


def isomers_c9h10n2o2pf2cl() -> GoalDirectedBenchmark:
    """
    Benchmark to try and get 100 isomers for C9H10N2O2PF2Cl.
    """

    specification = uniform_specification(100)

    return GoalDirectedBenchmark(name='C9H10N2O2PF2Cl',
                                 objective=IsomerScoringFunction('C9H10N2O2PF2Cl'),
                                 contribution_specification=specification)


def hard_cobimetinib() -> GoalDirectedBenchmark:
    smiles = 'OC1(CN(C1)C(=O)C1=C(NC2=C(F)C=C(I)C=C2)C(F)=C(F)C=C1)C1CCCCN1'

    modifier = ClippedScoreModifier(upper_x=0.7)
    os_tf = TanimotoScoringFunction(smiles, fp_type='FCFP4', score_modifier=modifier)
    os_ap = TanimotoScoringFunction(smiles, fp_type='ECFP6',
                                    score_modifier=MinGaussianModifier(mu=0.75, sigma=0.1))

    rot_b = RdkitScoringFunction(descriptor=num_rotatable_bonds,
                                 score_modifier=MinGaussianModifier(mu=3, sigma=1))

    rings = RdkitScoringFunction(descriptor=num_aromatic_rings,
                                 score_modifier=MaxGaussianModifier(mu=3, sigma=1))

    t_cns = ArithmeticMeanScoringFunction([os_tf, os_ap, rot_b, rings, CNS_MPO_ScoringFunction()])

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name='Cobimetinib MPO',
                                 objective=t_cns,
                                 contribution_specification=specification)


def hard_osimertinib() -> GoalDirectedBenchmark:
    smiles = 'COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc2nccc(n2)c3cn(C)c4ccccc34'

    modifier = ClippedScoreModifier(upper_x=0.8)
    similar_to_osimertinib = TanimotoScoringFunction(smiles, fp_type='FCFP4', score_modifier=modifier)

    but_not_too_similar = TanimotoScoringFunction(smiles, fp_type='ECFP6',
                                                  score_modifier=MinGaussianModifier(mu=0.85, sigma=0.1))

    tpsa_over_100 = RdkitScoringFunction(descriptor=tpsa,
                                         score_modifier=MaxGaussianModifier(mu=100, sigma=10))

    logP_scoring = RdkitScoringFunction(descriptor=logP,
                                        score_modifier=MinGaussianModifier(mu=1, sigma=1))

    make_osimertinib_great_again = ArithmeticMeanScoringFunction(
        [similar_to_osimertinib, but_not_too_similar, tpsa_over_100, logP_scoring])

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name='Osimertinib MPO',
                                 objective=make_osimertinib_great_again,
                                 contribution_specification=specification)


def hard_fexofenadine() -> GoalDirectedBenchmark:
    """
    make fexofenadine less greasy
    :return:
    """
    smiles = 'CC(C)(C(=O)O)c1ccc(cc1)C(O)CCCN2CCC(CC2)C(O)(c3ccccc3)c4ccccc4'

    modifier = ClippedScoreModifier(upper_x=0.8)
    similar_to_fexofenadine = TanimotoScoringFunction(smiles, fp_type='AP', score_modifier=modifier)

    tpsa_over_90 = RdkitScoringFunction(descriptor=tpsa,
                                        score_modifier=MaxGaussianModifier(mu=90, sigma=10))

    logP_under_4 = RdkitScoringFunction(descriptor=logP,
                                        score_modifier=MinGaussianModifier(mu=4, sigma=1))

    optimize_fexofenadine = ArithmeticMeanScoringFunction(
        [similar_to_fexofenadine, tpsa_over_90, logP_under_4])

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name='Fexofenadine MPO',
                                 objective=optimize_fexofenadine,
                                 contribution_specification=specification)


def start_pop_ranolazine() -> GoalDirectedBenchmark:
    ranolazine = 'COc1ccccc1OCC(O)CN2CCN(CC(=O)Nc3c(C)cccc3C)CC2'

    modifier = ClippedScoreModifier(upper_x=0.7)
    similar_to_ranolazine = TanimotoScoringFunction(ranolazine, fp_type='AP', score_modifier=modifier)

    logP_under_4 = RdkitScoringFunction(descriptor=logP,
                                        score_modifier=MaxGaussianModifier(mu=7, sigma=1))

    aroma = RdkitScoringFunction(descriptor=num_aromatic_rings,
                                 score_modifier=MinGaussianModifier(mu=1, sigma=1))

    fluorine = RdkitScoringFunction(descriptor=AtomCounter('F'),
                                    score_modifier=GaussianModifier(mu=1, sigma=1.0))

    optimize_ranolazine = ArithmeticMeanScoringFunction([similar_to_ranolazine, logP_under_4, fluorine, aroma])

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name='Ranolazine MPO',
                                 objective=optimize_ranolazine,
                                 contribution_specification=specification,
                                 starting_population=[ranolazine])


def weird_physchem() -> GoalDirectedBenchmark:
    min_bertz = RdkitScoringFunction(descriptor=bertz,
                                     score_modifier=MaxGaussianModifier(mu=1500, sigma=200))

    mol_under_400 = RdkitScoringFunction(descriptor=mol_weight,
                                         score_modifier=MinGaussianModifier(mu=400, sigma=40))

    aroma = RdkitScoringFunction(descriptor=num_aromatic_rings,
                                 score_modifier=MinGaussianModifier(mu=3, sigma=1))

    fluorine = RdkitScoringFunction(descriptor=AtomCounter('F'),
                                    score_modifier=GaussianModifier(mu=6, sigma=1.0))

    opt_weird = ArithmeticMeanScoringFunction(
        [min_bertz, mol_under_400, aroma, fluorine])

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name='Physchem MPO',
                                 objective=opt_weird,
                                 contribution_specification=specification)


def similarity_cns_mpo(smiles, molecule_name) -> GoalDirectedBenchmark:
    benchmark_name = f'{molecule_name}'
    os_tf = TanimotoScoringFunction(smiles, fp_type='FCFP4')
    os_ap = TanimotoScoringFunction(smiles, fp_type='AP')
    anti_fp = TanimotoScoringFunction(smiles, fp_type='ECFP6',
                                      score_modifier=MinGaussianModifier(mu=0.70, sigma=0.1))

    t_cns = ArithmeticMeanScoringFunction([os_tf, os_ap, anti_fp, CNS_MPO_ScoringFunction()])

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name=benchmark_name,
                                 objective=t_cns,
                                 contribution_specification=specification)


def similarity(smiles: str, name: str, fp_type: str = 'ECFP4', threshold: float = 0.7,
               rediscovery: bool = False) -> GoalDirectedBenchmark:
    category = 'rediscovery' if rediscovery else 'similarity'
    benchmark_name = f'{name} {category}'

    modifier = ClippedScoreModifier(upper_x=threshold)
    scoring_function = TanimotoScoringFunction(target=smiles, fp_type=fp_type, score_modifier=modifier)
    if rediscovery:
        specification = uniform_specification(1)
    else:
        specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name=benchmark_name,
                                 objective=scoring_function,
                                 contribution_specification=specification)


def logP_benchmark(target: float) -> GoalDirectedBenchmark:
    benchmark_name = f'logP (target: {target})'
    objective = RdkitScoringFunction(descriptor=logP,
                                     score_modifier=GaussianModifier(mu=target, sigma=1))

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name=benchmark_name,
                                 objective=objective,
                                 contribution_specification=specification)


def tpsa_benchmark(target: float) -> GoalDirectedBenchmark:
    benchmark_name = f'TPSA (target: {target})'
    objective = RdkitScoringFunction(descriptor=tpsa,
                                     score_modifier=GaussianModifier(mu=target, sigma=20.0))

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name=benchmark_name,
                                 objective=objective,
                                 contribution_specification=specification)


def cns_mpo() -> GoalDirectedBenchmark:
    specification = uniform_specification(1, 10, 100)
    return GoalDirectedBenchmark(name='CNS MPO', objective=CNS_MPO_ScoringFunction(),
                                 contribution_specification=specification)


def qed_benchmark() -> GoalDirectedBenchmark:
    specification = uniform_specification(1, 10, 100)
    return GoalDirectedBenchmark(name='QED',
                                 objective=RdkitScoringFunction(descriptor=qed),
                                 contribution_specification=specification)


def median_molecule() -> GoalDirectedBenchmark:
    t_camphor = TanimotoScoringFunction('CC1(C)C2CCC1(C)C(=O)C2', fp_type='ECFP4')
    t_menthol = TanimotoScoringFunction('CC(C)C1CCC(C)CC1O', fp_type='ECFP4')
    median = ArithmeticMeanScoringFunction([t_menthol, t_camphor])

    specification = uniform_specification(1, 10, 100)

    return GoalDirectedBenchmark(name='Median molecules',
                                 objective=median,
                                 contribution_specification=specification)


def novelty_benchmark(training_set_file: str, number_samples: int) -> DistributionLearningBenchmark:
    smiles_list = [s.strip() for s in open(training_set_file).readlines()]
    return NoveltyBenchmark(number_samples=number_samples, training_set=smiles_list)


def kldiv_benchmark(training_set_file: str, number_samples: int) -> DistributionLearningBenchmark:
    smiles_list = [s.strip() for s in open(training_set_file).readlines()]
    return KLDivBenchmark(number_samples=number_samples, training_set=smiles_list)


def frechet_benchmark(training_set_file: str, number_samples: int) -> DistributionLearningBenchmark:
    smiles_list = [s.strip() for s in open(training_set_file).readlines()]
    return FrechetBenchmark(training_set=smiles_list, sample_size=number_samples)
