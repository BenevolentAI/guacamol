import datetime
import json
import logging
from typing import List

from guacamol.distribution_learning_benchmark import DistributionLearningBenchmark
from guacamol.distribution_matching_generator import DistributionMatchingGenerator
from guacamol.benchmark_suites import distribution_learning_benchmark_suite

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def assess_distribution_learning(model: DistributionMatchingGenerator,
                                 chembl_training_file: str,
                                 json_output_file='output_distribution_learning.json',
                                 benchmark_version='v1') -> None:
    """
    Assesses a distribution-matching model for de novo molecule design.

    Args:
        model: Model to evaluate
        chembl_training_file: path to ChEMBL training set, necessary for some benchmarks
        json_output_file: Name of the file where to save the results in JSON format
        benchmark_version: which benchmark suite to execute
    """
    logger.info(f'Benchmarking distribution learning, version {benchmark_version}')
    benchmarks = distribution_learning_benchmark_suite(chembl_file_path=chembl_training_file,
                                                       version_name=benchmark_version)

    _evaluate_distribution_learning_benchmarks(model=model, benchmarks=benchmarks, json_output_file=json_output_file)


def _evaluate_distribution_learning_benchmarks(model: DistributionMatchingGenerator,
                                               benchmarks: List[DistributionLearningBenchmark],
                                               json_output_file: str) -> None:
    """
    Evaluate a model with the given benchmarks.
    Should not be called directly except for testing purposes.

    Args:
        model: model to assess
        benchmarks: list of benchmarks to evaluate
        json_output_file: Name of the file where to save the results in JSON format
    """

    logger.info(f'Number of benchmarks: {len(benchmarks)}')

    results = []
    for i, benchmark in enumerate(benchmarks, 1):
        logger.info(f'Running benchmark {i}/{len(benchmarks)}: {benchmark.name}')
        result = benchmark.assess_model(model)
        logger.info(f'Results for the benchmark "{result.benchmark_name}":')
        logger.info(f'  Score: {result.score:.6f}')
        logger.info(f'  Sampling time: {str(datetime.timedelta(seconds=int(result.sampling_time)))}')
        logger.info(f'  Metadata: {result.metadata}')
        results.append(result)

    logger.info('Finished execution of the benchmarks')

    all_results_dict = [vars(result) for result in results]

    logger.info(f'Save results to file {json_output_file}')
    with open(json_output_file, 'wt') as f:
        f.write(json.dumps(all_results_dict, indent=4, sort_keys=True))
