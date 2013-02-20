#/usr/bin/env python
import os.path as op
import argparse

from genelinkage.orfs_and_hmm import find_classify_orfs
from genelinkage.calc_tetra import generate_tetra
from genelinkage.calc_dists_and_plot import plot_dist

if __name__ == '__main__':

    desc = 'Find linkages between genes.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--working-dir', '-d', type=str, \
                        default='.',
                        help='Directory to save temp files in.')
    parser.add_argument('--project', '-p', type=str, \
                        default='',
                        help='Name of the project in the working directory.')
    parser.add_argument('--scaffolds', '-s', type=str, \
                        help='Path to the assembled contigs/scaffolds.')
    parser.add_argument('--markov-model', '-m', type=str, \
                        help='Path to the HMMER gene-prediction model.')

    parser.add_argument('--hmmer', type=str, default=None,
                        help='Location of HMMER 3.0')
    parser.add_argument('--metagenemark', type=str, default=None,
                        help='Location of MetaGeneMark')

    parser.add_argument('--all', '-a', action='store_true', \
                        dest='build_all', default=False, \
                        help='Run the entire toolchain.')
    parser.add_argument('--find-orfs', '-f', action='store_true', \
                        dest='build_orfs', default=False, \
                        help='Find and classify ORFs.')
    parser.add_argument('--tetra', '-t', action='store_true', \
                        dest='build_tetra', default=False, \
                        help='Calculate the tetranucleotide abundances.')
    parser.add_argument('--plot', action='store_true', \
                        dest='build_plot', default=False, \
                        help='Plot the gene-gene distances.')
    parser.add_argument('--dist', action='store_true', \
                        dest='build_dist', default=False, \
                        help='Save the gene-gene distances.')

    args = parser.parse_args()

    project = op.join(args.working_dir, args.project)
    if args.build_orfs or args.build_all:
        find_classify_orfs(project, args.scaffolds, args.markov_model, \
                           args.metagenemark, args.hmmer)

    if args.build_tetra or args.build_all:
        generate_tetra(project, args.scaffolds)

    if args.build_plot or args.build_dist or args.build_all:
        M_DRAWS = 10  # number of Monte Carlo draws to do for controls

        #gene_list = ['psaA', 'psaB', 'psbA', 'psbB', 'pufM', 'pufL', 'pr', 'pioA',
        #             'pioC', 'iro', 'coxB', 'ompC', 'arch_amoA', 'bact_amoA',
        #             'mmoZ', 'hszA', 'sqR-allo', 'sqR-rhodo', 'narG', 'nirK',
        #             'dsrA', 'dsrB', 'mcrA', 'frhB', 'cdhD', 'fdhA', 'mvK', 'dxr',
        #             'gggps', 'sqdB', 'cdsA-allo', 'cdsA-geo', 'cdsA-rhodo',
        #             'cdsA-synn', 'mglcD', 'mgdA', 'btaA', 'olsB', 'shc', 'osc',
        #             'cas1', 'crtI-allo', 'crtI-rhodo', 'crtP', 'nifH', 'luxI',
        #             'raiI', 'por', 'bchF', 'rpoB']
        gene_list = ['osc', 'shc', 'dsrA', 'dsrB']
        plot = args.build_plot or args.build_all
        save = args.build_dist or args.build_all
        plot_dist(project, gene_list, 30, M_DRAWS, plot, save)
