from guacamol.goal_directed_benchmark import GoalDirectedBenchmark
from guacamol.goal_directed_score_contributions import uniform_specification
from guacamol.similaritymetric_scoring_function import SimilarityMetricScoringFunction


def similarity_to_transition_state(method, ref_state, num_conformations,
                                   num_protonation_states, options):
    """
    Computes similarity score using a metric from the similaritymetrics repo.
    :param method:              string, name of the method to compute the metric
    :param ref_state:           string, name of the file with the reference transition state
    :param num_conformations:   int, number of conformations sampled
    :param options:             dict, options of the scorer.
    :return:                    GoalDirectedBenchmark
    """
    benchmark_name = "{} transition state similarity".format(ref_state.split(".")[0])
    scoring_function = SimilarityMetricScoringFunction(method, options, ref_state,
                                                       num_conformations,
                                                       num_protonation_states)
    specification = uniform_specification(1, 10, 100)
    return GoalDirectedBenchmark(name=benchmark_name,
                                 objective=scoring_function,
                                 contribution_specification=specification)

