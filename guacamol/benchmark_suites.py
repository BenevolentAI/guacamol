from typing import List

from sklearn.base import BaseEstimator

from guacamol.distribution_learning_benchmark import DistributionLearningBenchmark, ValidityBenchmark, UniquenessBenchmark
from guacamol.goal_directed_benchmark import GoalDirectedBenchmark
from guacamol.scoring_function import ArithmeticMeanScoringFunction
from guacamol.standard_benchmarks import hard_cobimetinib, similarity, logP_benchmark, cns_mpo, \
    qed_benchmark, median_camphor_menthol, novelty_benchmark, isomers_c11h24, isomers_c7h8n2o2, isomers_c9h10n2o2pf2cl, \
    frechet_benchmark, tpsa_benchmark, hard_osimertinib, hard_fexofenadine, weird_physchem, start_pop_ranolazine, \
    kldiv_benchmark, perindopril_rings, amlodipine_rings, sitagliptin_replacement, zaleplon_with_other_formula, valsartan_smarts, \
    median_tadalafil_sildenafil, decoration_hop, scaffold_hop, ranolazine_mpo, pioglitazone_mpo

def goal_directed_benchmark_suite(version_name: str) -> List[GoalDirectedBenchmark]:
    return GoalDirectedBenchmarkSuite().getBenchmark(version_name).getList()

# factory (as singleton software engineering pattern)
class GoalDirectedBenchmarkSuite(object):
    class __GoalDirectedBenchmarkSuite:
        def __init__(self):
            self.register = {}

        def getBenchmarkNames(self):
            return sorted(self.register.keys())
        
        def getBenchmark(self,benchmark_name ):
            if not len(benchmark_name):
                raise ValueError("getBenchmark called with no benchmark_name")
            benchmark=None
            try:
                benchmark = self.register[benchmark_name.lower()]
            except:
                error_msg="No Benchmark found named "+benchmark_name+"\nCurrently registered Benchmarks:\n\t"+"\n\t".join(self.register.keys())
                raise Exception(error_msg)
            return benchmark
    instance = None
    def __new__(cls): # __new__ always a classmethod
        if not GoalDirectedBenchmarkSuite.instance:
            GoalDirectedBenchmarkSuite.instance = GoalDirectedBenchmarkSuite.__GoalDirectedBenchmarkSuite()
        return GoalDirectedBenchmarkSuite.instance
    def __getattr__(self, name):
        return getattr(self.instance, name)
    def __setattr__(self, name):
        return setattr(self.instance, name)


# factory product
class AbstractFactoryProductGuacamol(BaseEstimator):
    def getHyperparameters(self):
        hyperparameters=self.get_params()
        return hyperparameters
    def getHyperparametersDescription(self):
        """Override me for the actual calculation"""
        raise NotImplementedError
    def getName(self):
        """Override me for the actual calculation"""
        raise NotImplementedError
    def getNameFromHyperparameters(self, prefix):
        description,abbreviation=self.getHyperparametersDescription()
        hyperparameters=self.getHyperparameters()
        
        _name=prefix
        for key in sorted (hyperparameters.keys()):
            _name+="-"
            _name+=abbreviation[key]
            _name+=str(hyperparameters[key])
        
        return _name.lower()


# factory product
class GoalDirectedBenchmarkSuiteEntry(AbstractFactoryProductGuacamol):
    def __init__(self):
        GoalDirectedBenchmarkSuite().register[self.getName()] = self

    def getList(self) -> List[GoalDirectedBenchmark]:
        """Override me for the actual benchmark list"""
        raise NotImplementedError


#define
class GoalDirectedV1(GoalDirectedBenchmarkSuiteEntry):
    def __init__(self, max_logP = 6.35584):
        self.max_logP=max_logP
        GoalDirectedBenchmarkSuiteEntry.__init__(self)

    def getName(self):
        _name='v1'
        return _name.lower()

    def getHyperparametersDescription(self):
        description={}
        abbreviation={}
        return description,abbreviation

    def getList(self) -> List[GoalDirectedBenchmark]:
        max_logP = self.max_logP
        return [
            isomers_c11h24(mean_function='arithmetic'),
            isomers_c7h8n2o2(mean_function='arithmetic'),
            isomers_c9h10n2o2pf2cl(mean_function='arithmetic', n_samples=100),
       
            hard_cobimetinib(max_logP=max_logP),
            hard_osimertinib(ArithmeticMeanScoringFunction),
            hard_fexofenadine(ArithmeticMeanScoringFunction),
            weird_physchem(),
       
            # start pop benchmark
            # e.g.
            start_pop_ranolazine(),
       
            # similarity Benchmarks
       
            # explicit rediscovery
            similarity(smiles='CC1=CC=C(C=C1)C1=CC(=NN1C1=CC=C(C=C1)S(N)(=O)=O)C(F)(F)F', name='Celecoxib', fp_type='ECFP4', threshold=1.0, rediscovery=True),
            similarity(smiles='Cc1c(C)c2OC(C)(COc3ccc(CC4SC(=O)NC4=O)cc3)CCc2c(C)c1O', name='Troglitazone', fp_type='ECFP4', threshold=1.0, rediscovery=True),
            similarity(smiles='CN(C)S(=O)(=O)c1ccc2Sc3ccccc3C(=CCCN4CCN(C)CC4)c2c1', name='Thiothixene', fp_type='ECFP4', threshold=1.0, rediscovery=True),
       
            # generate similar stuff
            similarity(smiles='Clc4cccc(N3CCN(CCCCOc2ccc1c(NC(=O)CC1)c2)CC3)c4Cl', name='Aripiprazole', fp_type='FCFP4', threshold=0.75),
            similarity(smiles='CC(C)(C)NCC(O)c1ccc(O)c(CO)c1', name='Albuterol', fp_type='FCFP4', threshold=0.75),
            similarity(smiles='COc1ccc2[C@H]3CC[C@@]4(C)[C@@H](CC[C@@]4(O)C#C)[C@@H]3CCc2c1', name='Mestranol', fp_type='AP', threshold=0.75),
       
            logP_benchmark(target=-1.0),
            logP_benchmark(target=8.0),
            tpsa_benchmark(target=150.0),
       
            cns_mpo(max_logP=max_logP),
            qed_benchmark(),
            median_camphor_menthol(ArithmeticMeanScoringFunction)
        ]

#register
GoalDirectedV1()

#define
class GoalDirectedV2(GoalDirectedBenchmarkSuiteEntry):
    def __init__(self):
        GoalDirectedBenchmarkSuiteEntry.__init__(self)

    def getName(self):
        _name='v2'
        return _name.lower()

    def getHyperparametersDescription(self):
        description={}
        abbreviation={}
        return description,abbreviation

    def getList(self) -> List[GoalDirectedBenchmark]:
        return [
            # explicit rediscovery
            similarity(smiles='CC1=CC=C(C=C1)C1=CC(=NN1C1=CC=C(C=C1)S(N)(=O)=O)C(F)(F)F', name='Celecoxib', fp_type='ECFP4', threshold=1.0, rediscovery=True),
            similarity(smiles='Cc1c(C)c2OC(C)(COc3ccc(CC4SC(=O)NC4=O)cc3)CCc2c(C)c1O', name='Troglitazone', fp_type='ECFP4', threshold=1.0, rediscovery=True),
            similarity(smiles='CN(C)S(=O)(=O)c1ccc2Sc3ccccc3C(=CCCN4CCN(C)CC4)c2c1', name='Thiothixene', fp_type='ECFP4', threshold=1.0, rediscovery=True),
       
            # generate similar stuff
            similarity(smiles='Clc4cccc(N3CCN(CCCCOc2ccc1c(NC(=O)CC1)c2)CC3)c4Cl', name='Aripiprazole', fp_type='ECFP4', threshold=0.75),
            similarity(smiles='CC(C)(C)NCC(O)c1ccc(O)c(CO)c1', name='Albuterol', fp_type='FCFP4', threshold=0.75),
            similarity(smiles='COc1ccc2[C@H]3CC[C@@]4(C)[C@@H](CC[C@@]4(O)C#C)[C@@H]3CCc2c1', name='Mestranol', fp_type='AP', threshold=0.75),
       
            # isomers
            isomers_c11h24(),
            isomers_c9h10n2o2pf2cl(),
       
            # median molecules
            median_camphor_menthol(),
            median_tadalafil_sildenafil(),
       
            # all other MPOs
            hard_osimertinib(),
            hard_fexofenadine(),
            ranolazine_mpo(),
            perindopril_rings(),
            amlodipine_rings(),
            sitagliptin_replacement(),
            zaleplon_with_other_formula(),
            valsartan_smarts(),
            decoration_hop(),
            scaffold_hop(),
        ]
    
#register
GoalDirectedV2()

#define
class GoalDirectedTrivial(GoalDirectedBenchmarkSuiteEntry):
    def __init__(self):
        GoalDirectedBenchmarkSuiteEntry.__init__(self)

    def getName(self):
        _name='trivial'
        return _name.lower()

    def getHyperparametersDescription(self):
        description={}
        abbreviation={}
        return description,abbreviation

    def getList(self) -> List[GoalDirectedBenchmark]:
        """
        Trivial goal-directed benchmarks from the paper.
        """
        return [
            logP_benchmark(target=-1.0),
            logP_benchmark(target=8.0),
            tpsa_benchmark(target=150.0),
            cns_mpo(),
            qed_benchmark(),
            isomers_c7h8n2o2(),
            pioglitazone_mpo(),
        ]

#register
GoalDirectedTrivial()


#define
class GoalDirectedAll(GoalDirectedBenchmarkSuiteEntry):
    def __init__(self, max_logP = 6.35584):
        self.max_logP=max_logP
        GoalDirectedBenchmarkSuiteEntry.__init__(self)

    def getName(self):
        _name='all'
        return _name.lower()

    def getHyperparametersDescription(self):
        description={}
        abbreviation={}
        return description,abbreviation

    def getList(self) -> List[GoalDirectedBenchmark]:
        max_logP = self.max_logP
        return [
            # explicit rediscovery
            similarity(smiles='CC1=CC=C(C=C1)C1=CC(=NN1C1=CC=C(C=C1)S(N)(=O)=O)C(F)(F)F', name='Celecoxib', fp_type='ECFP4', threshold=1.0, rediscovery=True),
            similarity(smiles='Cc1c(C)c2OC(C)(COc3ccc(CC4SC(=O)NC4=O)cc3)CCc2c(C)c1O', name='Troglitazone', fp_type='ECFP4', threshold=1.0, rediscovery=True),
            similarity(smiles='CN(C)S(=O)(=O)c1ccc2Sc3ccccc3C(=CCCN4CCN(C)CC4)c2c1', name='Thiothixene', fp_type='ECFP4', threshold=1.0, rediscovery=True),
       
            # generate similar stuff
            similarity(smiles='Clc4cccc(N3CCN(CCCCOc2ccc1c(NC(=O)CC1)c2)CC3)c4Cl', name='Aripiprazole', fp_type='ECFP4', threshold=0.75),
            similarity(smiles='CC(C)(C)NCC(O)c1ccc(O)c(CO)c1', name='Albuterol', fp_type='FCFP4', threshold=0.75),
            similarity(smiles='COc1ccc2[C@H]3CC[C@@]4(C)[C@@H](CC[C@@]4(O)C#C)[C@@H]3CCc2c1', name='Mestranol', fp_type='AP', threshold=0.75),
       
            # isomers
            isomers_c7h8n2o2(),
            isomers_c7h8n2o2(mean_function='arithmetic'),
            isomers_c11h24(),
            isomers_c11h24(mean_function='arithmetic'),
            isomers_c9h10n2o2pf2cl(),
            isomers_c9h10n2o2pf2cl(mean_function='arithmetic', n_samples=100),
       
            # median molecules
            median_camphor_menthol(),
            median_camphor_menthol(ArithmeticMeanScoringFunction),
            median_tadalafil_sildenafil(),
       
            # start pop benchmark
            # e.g.
            start_pop_ranolazine(),

            # scaffold and decoration hopping
            decoration_hop(),
            scaffold_hop(),

            # all other MPOs
            amlodipine_rings(),
            cns_mpo(),
            cns_mpo(max_logP=max_logP),
            hard_cobimetinib(max_logP=max_logP),
            hard_fexofenadine(),
            hard_fexofenadine(ArithmeticMeanScoringFunction),
            hard_osimertinib(),
            hard_osimertinib(ArithmeticMeanScoringFunction),
            logP_benchmark(target=-1.0),
            logP_benchmark(target=8.0),
            perindopril_rings(),
            pioglitazone_mpo(),
            qed_benchmark(),
            ranolazine_mpo(),
            sitagliptin_replacement(),
            tpsa_benchmark(target=150.0),
            valsartan_smarts(),
            weird_physchem(),
            zaleplon_with_other_formula(),
        ]
    
#register
GoalDirectedAll()


# factory (as singleton software engineering pattern)
class DistributionLearningBenchmarkSuite(object):
    class __DistributionLearningBenchmarkSuite:
        def __init__(self):
            self.register = {}

        def getBenchmarkNames(self):
            return sorted(self.register.keys())
        
        def getBenchmark(self,benchmark_name ):
            if not len(benchmark_name):
                raise ValueError("getBenchmark called with no benchmark_name")
            benchmark=None
            try:
                benchmark = self.register[benchmark_name.lower()]
            except:
                error_msg="No Benchmark found named "+benchmark_name+"\nCurrently registered Benchmarks:\n\t"+"\n\t".join(self.register.keys())
                raise Exception(error_msg)
            return benchmark
    instance = None
    def __new__(cls): # __new__ always a classmethod
        if not DistributionLearningBenchmarkSuite.instance:
            DistributionLearningBenchmarkSuite.instance = DistributionLearningBenchmarkSuite.__DistributionLearningBenchmarkSuite()
        return DistributionLearningBenchmarkSuite.instance
    def __getattr__(self, name):
        return getattr(self.instance, name)
    def __setattr__(self, name):
        return setattr(self.instance, name)


def distribution_learning_benchmark_suite(chembl_file_path: str,
                                          version_name: str,
                                          number_samples: int) -> List[DistributionLearningBenchmark]:
    """
    Returns a suite of benchmarks for a specified benchmark version

    Args:
        chembl_file_path: path to ChEMBL training set, necessary for some benchmarks
        version_name: benchmark version

    Returns:
        List of benchmaks
    """

    # For distribution-learning, v1 and v2 are identical
    #register, of not already registered
    if not (version_name in DistributionLearningBenchmarkSuite().getBenchmarkNames()):
        DistributionLearningV1(chembl_file_path=chembl_file_path, number_samples=number_samples)
        DistributionLearningV1(registerAsName='v2',chembl_file_path=chembl_file_path, number_samples=number_samples)

    #return
    return DistributionLearningBenchmarkSuite().getBenchmark(version_name).getList()

# factory product
class DistributionLearningBenchmarkSuiteEntry(AbstractFactoryProductGuacamol):
    def __init__(self):
        DistributionLearningBenchmarkSuite().register[self.getName()] = self

    def getList(self) -> List[DistributionLearningBenchmark]:
        """Override me for the actual benchmark list"""
        raise NotImplementedError


#define
class DistributionLearningV1(DistributionLearningBenchmarkSuiteEntry):
    def __init__(self, registerAsName=None, chembl_file_path: str = None, number_samples: int = 10000, shortname="chembl"):
        self.registerAsName=registerAsName
        self.chembl_file_path=chembl_file_path
        self.number_samples=number_samples
        self.shortname=shortname
        DistributionLearningBenchmarkSuiteEntry.__init__(self)

    def getName(self):
        _name='v1'
        if self.registerAsName is not None:
            _name=self.registerAsName
        if self.shortname is not None:
            _name+='-'
            _name+=self.shortname
        return _name.lower()

    def getHyperparametersDescription(self):
        description={}
        abbreviation={}
        return description,abbreviation

    def getList(self) -> List[DistributionLearningBenchmark]:
        """
        Suite of distribution learning benchmarks, v1.
    
        Args:
            chembl_file_path: path to the file with the reference ChEMBL molecules
    
        Returns:
            List of benchmarks, version 1
        """
        return [
            ValidityBenchmark(number_samples=self.number_samples),
            UniquenessBenchmark(number_samples=self.number_samples),
            novelty_benchmark(training_set_file=self.chembl_file_path, number_samples=self.number_samples),
            kldiv_benchmark(training_set_file=self.chembl_file_path, number_samples=self.number_samples),
            frechet_benchmark(training_set_file=self.chembl_file_path, number_samples=self.number_samples)
        ]
    
