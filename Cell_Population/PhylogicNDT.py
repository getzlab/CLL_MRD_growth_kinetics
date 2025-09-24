import argparse
import logging
import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/")

# Remove all handlers associated with the root logger object.
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'phylogicndt.log')
print(filename)
logging.basicConfig(filename=filename,
                    filemode='w',
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',
                    level=getattr(logging, "INFO"))

import CellPopulation

def build_parser():
    parser = argparse.ArgumentParser(description="Run PhylogicNDT")

    # global parameters common to most tools
    base_parser = argparse.ArgumentParser(add_help=False)

    # option for specifying individual/patient ID
    base_parser.add_argument('--indiv_id', '-i',
                             type=str,
                             action='store',
                             dest='indiv_id',
                             default='Indiv1',
                             help='Patient/Case ID')

    # Samples information
    # Specifying samples on cmdline one-by-one in format sample_id:maf_fn:seg_fn:purity:timepoint
    base_parser.add_argument("-s", "--sample", dest='sample_data', action='append', type=str,
                             help="Sample data, format sample_id:maf_fn:seg_fn:purity:timepoint; each sample specify separately by multiple -s ..; ")

    # Alternative: Instead of specifying samples on cmdline, a tsv - sif (sample information file) may be used.
    base_parser.add_argument("-sif", "--sample_information_file", dest='sif', type=str, help="""Sample information tsv file with sample_ids and CCF and copy-number file_paths; \n 
    format per row (with header) sample_id\tmaf_fn\tseg_fn\tpurity\ttimepoint""")

    # Filtering of Mutations/Events
    # option for specifying blacklist
    base_parser.add_argument("-bl", '--artifact_blacklist',
                             type=str,
                             action='store',
                             dest='artifact_blacklist',
                             default=os.path.join(os.path.dirname(__file__), 'data/supplement_data/Blacklist_SNVs.txt'),
                             help='path to blacklist')
    # whitelist - always overrules blacklist.
    base_parser.add_argument("-wl", '--artifact_whitelist',
                             type=str,
                             action='store',
                             dest='artifact_whitelist',
                             default='',
                             help='path to artifact whitelist - takes precedence over blacklist')

    # option for specifying custom drivers list
    base_parser.add_argument("-drv", '--driver_genes_file',
                             type=str,
                             action='store',
                             dest='driver_genes_file',
                             default=os.path.join(os.path.dirname(__file__),
                                                  'data/supplement_data/Driver_genes_v1.0.txt'),
                             help='driver list file')

    base_parser.add_argument("-tr", '--treatment_data',
                             type=str,
                             action='store',
                             dest='treatment_data',
                             default=None,
                             help='path to treatment data file')

    base_parser.add_argument('-ts', '--tumor_size',
                             action='store',
                             dest='tumor_size',
                             default=None)

    base_parser.add_argument('--blacklist_threshold', '-bt',
                             type=float,
                             action='store',
                             dest='blacklist_threshold',
                             default=0.1,
                             help='ccf threshold for blacklisting clusters for a BuildTree and Cell Population')

    base_parser.add_argument('--seed',
                             type=int,
                             action='store',
                             dest='seed',
                             default=None,
                             help='input a random seed for reproducibility')

    # different Tools of the PhylogicNDT Package
    subparsers = parser.add_subparsers(title="tool", description="Choose a tool to run", help='Try the Cluster tool')

    cellpopulation = subparsers.add_parser("CellPopulation", help="CellPopulation module for computing cell abundance.",
                                           parents=[base_parser])
    cellpopulation.add_argument('--cluster_ccf', '-c',
                                type=str,
                                action='store',
                                dest='cluster_ccf_file',
                                help='tsv file phylogic clustering results')
    cellpopulation.add_argument('--mutation_ccf', '-m',
                                type=str,
                                action='store',
                                dest='mutation_ccf_file',
                                help='tsv file generated by clustering')
    cellpopulation.add_argument('--tree_tsv', '-t',
                                type=str,
                                action='store',
                                dest='tree_tsv',
                                help='tsv file generated by build tree module')
    cellpopulation.add_argument('--n_iter', '-ni',
                                type=int,
                                action='store',
                                dest='n_iter',
                                default=250,
                                help='number iterations')

    cellpopulation.add_argument('--tree_number',
                                dest='tree_number',
                                action='store',
                                type=int,
                                help='Specify which tree to select from the ranked build_tree_posteriors (1-indexed); default is 1 (most likely tree).',
                                default=1)
    cellpopulation.set_defaults(func=CellPopulation.run_tool)

    return parser.parse_args()

# if __name__ == "main":
parser = build_parser()
args = parser
args.func(args)