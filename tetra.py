#/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import os
import shutil
import itertools
import subprocess
#from collections import Counter
from glob import iglob
from functools import partial
from multiprocessing import Pool

import numpy as np
from scipy.stats import gaussian_kde, mannwhitneyu
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

DATA_DIR = "/n/pearsonfs1/SequenceData/MahoneyLake7M/"
GENE_DIR = '/n/pearsonfs1/Roderick/create_COGs/Genes'
GFF_FILE = DATA_DIR + "Genes/metagenemark.gff"
GFA_FILE = DATA_DIR + "Genes/genes.faa"
FASTA_FILE = DATA_DIR + 'scaffold.fa'
ALL_GENES = GENE_DIR + "all.sto"
GENE_DB = GENE_DIR + "genes.hmm"
GENE_FILE = "../MHL7_hmm_results.txt"
TETRA_FILE = "../MHL7_hmm_tetra.txt"
DIST_FILE = "../MHL7_hmm_dist.txt"

build_orfs = False
build_hmmdb = False
build_tetra = False
build_dist = True


if build_orfs:
    #TODO: META_MOD should point to actual MetaGeneMark DB
    META_MOD = "./MetaGeneMark_linux64/MetaGeneMark_v1.mod"
    subprocess.call(['gmhmmp', '-a', '-f', 'G', '-m', META_MOD, '-o',
                     GFF_FILE, FASTA_FILE])
    subprocess.call(['aa_from_gff.pl', '<', GFF_FILE, '>', GFA_FILE])

if build_hmmdb:
    #TODO: do the blasting on all the individual gene *.FA files
    # to assemble the aligned sequence DBs for each one
    all_gene_file = open(ALL_GENES, 'w')
    for filename in iglob(os.path.join(GENE_DIR, '*.sto')):
        with open(filename, 'r') as f:
            shutil.copyfileobj(f, ALL_GENES)
    all_gene_file.close()

    subprocess.call(["hmmbuild", GENE_DB, ALL_GENES])
    subprocess.call(["hmmpress", GENE_DB])
    subprocess.call(["hmmscan", "--tblout", GENE_FILE, GENE_DB, GFA_FILE])

if build_tetra:
    import screed

    def rc(seq):  # reverse complement
        """
        Function finds reverse complement of a sequence.
        """
        invert_table = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(invert_table.get(i, 'N') for i in seq[::-1])

    def slid_win(seq, size=4):
        """Returns a sliding window along seq."""
        itr = iter(seq)
        buf = ''
        for _ in range(size):
            buf += next(itr)
        for l in itr:
            yield buf
            buf = buf[1:] + l
        yield buf

    # generate a mapping to merge reverse complement tetranucleotides
    seq_map = {'': 4 * 'N'}
    for s in (''.join(i) for i in itertools.product(*(4 * ['ATGC']))):
        if rc(s) not in seq_map or s not in seq_map:
            seq_map[s] = s
            seq_map[rc(s)] = s

    # write out the headers for the gene/tetranucleotide file
    fout = open(TETRA_FILE, 'w')
    srted_vals = list(set(seq_map.values()))
    srted_vals.sort()
    fout.write(','.join(['Gene'] + srted_vals) + '\n')

    f = open(GFF_FILE, "r")
    for i in range(7):
        f.readline()

    # make a dictionary to map genes back to their contigs
    gene2contig = {}
    for ln in f:
        flds = ln.strip().split('\t')
        gene2contig[flds[-1].replace(' ', '_')] = flds[0].split(' ')[0]
    f.close()

    hmm = open(GENE_FILE, 'r')
    contigs = screed.read_fasta_sequences(FASTA_FILE)

    for i in range(3):
        hmm.readline()

    p_gff_name = ''
    for ln in hmm:
        gene_name = ln[0:21].strip()
        gff_name = ln[32:53].strip()

        # if the same gene_id is listed multiple times in a row
        # that means HMMER found multiple matches for it. We only
        # want the first one (with the lowest E-value).
        if gff_name != p_gff_name:
            p_gff_name = gff_name
            frq = dict([(s, 0) for s in seq_map.values()])
            cc = contigs[gene2contig[gff_name]]
            for ss in slid_win(str(cc.sequence).upper(), 4):
                frq[seq_map.get(ss, 'NNNN')] += 1

            sum_frq = float(sum(frq.values()))
            if sum_frq == 0:
                sum_frq = 1
            fout.write(','.join([gene_name] + [str(frq[i] / sum_frq) \
                                            for i in srted_vals]))
            fout.write('\n')
            fout.flush()
    hmm.close()
    fout.close()

# start computing gene distances
if build_dist:
    M_DRAWS = 1000  # number of Monte Carlo draws to do
    G_SAMP = 1000  # number of each of both genes to sample in each draw
    FIG_FILE = '../MHL7_tetra.png'

    # convenience functions to return a gene name or
    # its tetranucleotide array
    up_name = lambda l: l.split(',')[0].strip()
    up_tet = lambda l: np.array([float(i) for i in l.split(',')[1:]])

    # first, get a count of how many of each gene there are
    with open(TETRA_FILE, 'r') as f:
        # how many dimensions are in each tetra score?
        tetlen = len(f.readline().split(',')) - 1
        gene_names, gene_ids = {}, []
        gene_ct = Counter()
        for l in f:
            gene_name = up_name(l)
            gene_ct[gene_name] += 1
            if gene_name not in gene_names:
                gene_names[gene_name] = len(gene_names)
            gene_ids.append(gene_names[gene_name])
        gene_ids = np.array(gene_ids)
        f.seek(0)
        f.readline()
        gene_tetra = np.empty((len(gene_ids), tetlen))
        for i, l in enumerate(f):
            gene_tetra[i] = up_tet(l)
    # get a count of how many of each gene there are
    print('Done buffering')

    def rnd_genes(genes=[], n=1):
        """
        Pick n random HMMER-identified "genes" with replacement
        from the gene-tetranucleotide file.

        If genes == [] pick from all genes.
        """
        # how many genes are there total?
        if genes == []:
            sel_genes = np.ones(gene_ids.shape, dtype=bool)
        else:
            sel_genes = np.zeros(gene_ids.shape, dtype=bool)
            for gene in genes:
                sel_genes = np.logical_or(sel_genes, \
                                          gene_ids == gene_names[gene])
        # randomly pick genes from the collection
        rand_picks = np.random.randint(sum(sel_genes), size=(n,))
        tetra = gene_tetra[sel_genes][rand_picks]
        return tetra

    def min_dst(tet1, tet2, allow_zero=True):
        dists = np.empty(tet1.shape[0])
        for i, t1 in enumerate(tet1):
            min_dist = np.sum((tet2 - t1) ** 2, axis=1)
            if not allow_zero:
                dists[i] = np.min(min_dist[min_dist != 0])
            else:
                dists[i] = np.min(min_dist)
        return dists

    def samp_dist(gene2, gene1):
        tet1 = rnd_genes([gene1], G_SAMP)
        if gene2 is None:
            g2 = []
        else:
            g2 = [gene2]
        tet2 = rnd_genes(g2, G_SAMP)
        dist = []
        for _ in range(M_DRAWS):
            dist += list(min_dst(tet1, tet2, allow_zero=(gene1 != gene2)))
        return dist

    #TODO: no cerI in this list?
    gene_list = ['psaA', 'psaB', 'psbA', 'psbB', 'pufM', 'pufL', 'pr', 'pioA',
                 'pioC', 'iro', 'coxB', 'ompC', 'arch_amoA', 'bact_amoA',
                 'mmoZ', 'hszA', 'sqR-allo', 'sqR-rhodo', 'narG', 'nirK',
                 'dsrA', 'dsrB', 'mcrA', 'frhB', 'cdhD', 'fdhA', 'mvK', 'dxr',
                 'gggps', 'sqdB', 'cdsA-allo', 'cdsA-geo', 'cdsA-rhodo',
                 'cdsA-synn', 'mglcD', 'mgdA', 'btaA', 'olsB', 'shc', 'osc',
                 'cas1', 'crtI-allo', 'crtI-rhodo', 'crtP', 'nifH', 'luxI',
                 'raiI', 'por', 'bchF', 'rpoA', 'rpoB']
    gene_list = [g for g in gene_list if gene_ct[g] >= 30]
    #gene_list = ['osc', 'shc', 'dsrA', 'dsrB']
    #gene_list = all_genes

    po = Pool()

    gs = gridspec.GridSpec(len(gene_list), len(gene_list))
    gs.update(wspace=0, hspace=0)
    xs = np.linspace(0, 0.005, 200)

    #dist_f = partial(samp_dist, gene1='osc')
    #dists = po.map(dist_f, ['shc'])

    for i, g1 in enumerate(gene_list):
        dist_f = partial(samp_dist, gene1=g1)
        dists = po.map(dist_f, gene_list)
        ctrls = po.map(dist_f, [None] * len(gene_list))
        for j, g2 in enumerate(gene_list):
            ax = plt.subplot(gs[i + j * len(gene_list)])
            mwu = mannwhitneyu(ctrls[j], dists[j])[1]
            txt = g1 + '->' + g2 + '\np={:.2e}'.format(mwu)
            ax.text(0.5, 0.95, txt, fontsize=2, \
              va='top', ha='center', transform=ax.transAxes)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            if sum(dists[j]) != 0:
                #print(g1, g2, dists[j])
                ax.plot(xs, gaussian_kde(dists[j])(xs), 'k-')
            if sum(ctrls[j]) != 0:
                ax.plot(xs, gaussian_kde(ctrls[j])(xs), 'r-')
            #TODO: scipy.stats.mannwhitneyu(dist, ctrl)
    plt.gcf().set_size_inches(24, 24)
    plt.savefig(FIG_FILE, dpi=300, bbox_inches='tight')
    #plt.show()

# library(ggplot2)
# d <- read.csv('gene_dist.txt')
# genes <- c('dsrA','dsrB','dsrC','dsrL','dsrN')
# ggplot(subset(d, FromGene %in% c(genes) & ToGene %in% c(genes)),
#   aes(x=Dist)) + facet_grid(FromGene ~ ToGene, labeller=label_both)
#   + geom_density() + xlim(0,1)
