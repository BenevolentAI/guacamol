import argparse
import gzip
import hashlib
import os.path
import pkgutil
import platform
from typing import List, Iterable

import numpy as np
from joblib import Parallel, delayed

from guacamol.utils.chemistry import canonicalize_list, filter_and_canonicalize, \
    initialise_neutralisation_reactions, split_charged_mol, get_fingerprints_from_smileslist
from guacamol.utils.data import download_if_not_present, get_time_string

TRAIN_HASH = '05ad85d871958a05c02ab51a4fde8530'
VALID_HASH = 'e53db4bff7dc4784123ae6df72e3b1f0'
TEST_HASH = '677b757ccec4809febd83850b43e1616'

CHEMBL_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_24_1/chembl_24_1_chemreps.txt.gz'
CHEMBL_FILE_NAME = 'chembl_24_1_chemreps.txt.gz'


def get_argparser():
    timestring = get_time_string()
    parser = argparse.ArgumentParser(description='Data Preparation for GuacaMol',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--destination', default='.', help='Download and Output location')
    parser.add_argument('-i', '--input', default=None, help='Filename of input smiles file')
    parser.add_argument('--output_prefix', default=timestring, help='Prefix of the output file')
    parser.add_argument('--n_jobs', default=8, type=int, help='Number of cores to use')
    parser.add_argument('--tanimoto_cutoff', default=0.323, type=float,
                        help='Remove molecules too similar to the holdout set')
    parser.add_argument('--chembl', action='store_true',
                        help='Specify to download and process molecules from chembl')
    return parser


def extract_chembl(line) -> str:
    """
    Extract smiles from chembl tsv

    Returns:
        SMILES string
    """
    return line.split('\t')[1]


def extract_smilesfile(line) -> str:
    """
    Extract smiles from SMILES file

    Returns:
        SMILES string
    """
    return line.split(' ')[0].strip()


class AllowedSmilesCharDictionary(object):
    """
    A fixed dictionary for druglike SMILES.
    """

    def __init__(self) -> None:
        self.forbidden_symbols = {'Ag', 'Al', 'Am', 'Ar', 'At', 'Au', 'D', 'E', 'Fe', 'G', 'K', 'L', 'M', 'Ra', 'Re',
                                  'Rf', 'Rg', 'Rh', 'Ru', 'T', 'U', 'V', 'W', 'Xe',
                                  'Y', 'Zr', 'a', 'd', 'f', 'g', 'h', 'k', 'm', 'si', 't', 'te', 'u', 'v', 'y'}

    def allowed(self, smiles: str) -> bool:
        """
        Determine if SMILES string has illegal symbols

        Args:
            smiles: SMILES string

        Returns:
            True if all legal
        """
        for symbol in self.forbidden_symbols:
            if symbol in smiles:
                print('Forbidden symbol {:<2}  in  {}'.format(symbol, smiles))
                return False
        return True


def get_raw_smiles(file_name, smiles_char_dict, open_fn, extract_fn) -> List[str]:
    """
    Extracts the raw smiles from an input file.
    open_fn will open the file to iterate over it (e.g. use open_fn=open or open_fn=filegzip.open)
    extract_fn specifies how to process the lines, choose from
    Pre-filter molecules of 5 <= length <= 200, because processing larger molecules (e.g. peptides) takes very long.

    Returns:
       a list of SMILES strings
    """
    data = []
    # open the gzipped chembl filegzip.open
    with open_fn(file_name, 'rt') as f:

        line_count = 0
        for line in f:

            line_count += 1
            # extract the canonical smiles column
            if platform.system() == "Windows":
                line = line.decode("utf-8")

            # smiles = line.split('\t')[1]

            smiles = extract_fn(line)

            # only keep reasonably sized molecules
            if 5 <= len(smiles) <= 200:

                smiles = split_charged_mol(smiles)

                if smiles_char_dict.allowed(smiles):
                    # check whether the molecular graph consists of
                    # multiple connected components (eg. in salts)
                    # if so, just keep the largest one

                    data.append(smiles)

        print(f'Processed {line_count} lines.')

    return data


def write_smiles(dataset: Iterable[str], filename: str):
    """
    Dumps a list of SMILES into a file, one per line
    """
    n_lines = 0
    with open(filename, 'w') as out:
        for smiles_str in dataset:
            out.write('%s\n' % smiles_str)
            n_lines += 1
    print(f'{filename} contains {n_lines} molecules')


def compare_hash(output_file: str, correct_hash: str):
    """
    Computes the md5 hash of a SMILES file and check it against a given one
    Raises an exception if hashes are different
    """
    output_hash = hashlib.md5(open(output_file, 'rb').read()).hexdigest()
    if output_hash != correct_hash:
        # print(f'{output_file} file has different hash {output_hash} than expected {correct_hash}!')
        raise ValueError(f'{output_file} file has different hash {output_hash} than expected {correct_hash}!')


def main():
    """ Get Chembl-23.

    Preprocessing steps:

    1) filter SMILES shorter than 5 and longer than 200 chars and those with forbidden symbols
    2) canonicalize, neutralize, only permit smiles shorter than 100 chars
    3) shuffle, write files, check if they are consistently hashed.
    """
    argparser = get_argparser()
    args = argparser.parse_args()

    # Set constants
    np.random.seed(1337)
    neutralization_rxns = initialise_neutralisation_reactions()
    smiles_dict = AllowedSmilesCharDictionary()

    tanimoto_cutoff = args.tanimoto_cutoff

    # Either use chembl, or supplied SMILES file.

    print('Preprocessing molecules...')

    if args.chembl:

        print('Using Chembl')

        chembl_file = os.path.join(args.destination, CHEMBL_FILE_NAME)

        data = pkgutil.get_data('guacamol.data', 'holdout_set_gcm_v1.smiles').decode('utf-8').splitlines()

        holdout_mols = [i.split(' ')[0] for i in data]
        holdout_set = set(canonicalize_list(holdout_mols, False))
        holdout_fps = get_fingerprints_from_smileslist(holdout_set)

        # Download Chembl23 if needed.
        download_if_not_present(chembl_file,
                                uri=CHEMBL_URL)
        raw_smiles = get_raw_smiles(chembl_file, smiles_char_dict=smiles_dict, open_fn=gzip.open,
                                    extract_fn=extract_chembl)

        file_prefix = 'chembl24_canon'

        print(f'Excluding molecules based on ECFP4 similarity of > {tanimoto_cutoff} to the holdout set')

    else:
        if args.input is None:
            raise IOError(
                'You need to specify an input smiles file with -i {file} or --input {file}. \n'
                'Alternatively, provide the --chembl flag to download and process molecules from ChEMBL24 (recommended)')

        raw_smiles = get_raw_smiles(args.input, smiles_char_dict=smiles_dict, open_fn=open,
                                    extract_fn=extract_smilesfile)
        tanimoto_cutoff = 100  # effectively no cutoff
        holdout_set = set([])
        holdout_fps = []
        file_prefix = args.output_prefix

    print()
    print(f'Standardizing {len(raw_smiles)} molecules using {args.n_jobs} cores...')

    # Process all the SMILES in parallel
    runner = Parallel(n_jobs=args.n_jobs, verbose=2)

    joblist = (delayed(filter_and_canonicalize)(smiles_str,
                                                holdout_set,
                                                holdout_fps,
                                                neutralization_rxns,
                                                tanimoto_cutoff,
                                                False)
               for smiles_str in raw_smiles)

    output = runner(joblist)

    # Put all nonzero molecules in a list, remove duplicates, sort and shuffle

    all_good_mols = sorted(list(set([item[0] for item in output if item])))
    np.random.shuffle(all_good_mols)
    print(f'Ended up with {len(all_good_mols)} molecules. Preparing splits...')

    # Split into train-dev-test
    # Check whether the md5-hashes of the generated smiles files match
    # the precomputed hashes, this ensures everyone works with the same splits.

    VALID_SIZE = int(0.05 * len(all_good_mols))
    TEST_SIZE = int(0.15 * len(all_good_mols))

    dev_set = all_good_mols[0:VALID_SIZE]
    dev_path = os.path.join(args.destination, f'{file_prefix}_dev-valid.smiles')
    write_smiles(dev_set, dev_path)

    test_set = all_good_mols[VALID_SIZE:VALID_SIZE + TEST_SIZE]
    test_path = os.path.join(args.destination, f'{file_prefix}_test.smiles')
    write_smiles(test_set, test_path)

    train_set = all_good_mols[VALID_SIZE + TEST_SIZE:]
    train_path = os.path.join(args.destination, f'{file_prefix}_train.smiles')
    write_smiles(train_set, train_path)

    # for chembl, check the hashes
    if args.chembl:
        compare_hash(train_path, TRAIN_HASH)
        compare_hash(dev_path, VALID_HASH)
        compare_hash(test_path, TEST_HASH)

        print('The train/test/dev-file md5 hashes match the expected hashes.')

    print('You are ready to go.')


if __name__ == '__main__':
    main()
