#/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import os
import shutil
import itertools
import subprocess
from collections import Counter
from glob import iglob

import screed
import numpy as np
from scipy.stats import gaussian_kde
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
build_tetra = True
build_dist = False


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
    def rc(seq):  # reverse complement
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
    FIG_FILE = './MHL7_tetra.png'

    up_name = lambda l: l.split(',')[0].strip()
    up_tet = lambda l: [float(i) for i in l.split(',')[1:]]

    with open(TETRA_FILE, 'r') as f:
        # how many dimensions are in each tetra score?
        tetlen = len(f.readline().split(',')) - 1

        # get a count of how many of each gene there are
        gene_ct = Counter(up_name(l) for l in f)
    all_genes = set(gene_ct)

    def rand_genes(genes, n=1):
        # how many genes are there total?
        ngenes = sum(gene_ct[g] for g in set(genes))
        # randomly pick genes from the collection
        gene_nums = np.random.randint(ngenes, size=(n,))
        with open(TETRA_FILE, 'r') as f:
            f.readline()  # skip the header line
            c = 0
            tetra = np.empty((n, tetlen), dtype=float)
            for l in f:
                if up_name(l) in genes:
                    tetra[np.where(gene_nums == c)[0]] = up_tet(l)
                    c += 1
        return tetra

    def min_dist(tet, genes, allow_zero=True):
        with open(TETRA_FILE, 'r') as f:
            f.readline()  # skip the header line
            dist = np.inf * np.ones(len(tet))
            for l in f:
                if up_name(l) in genes:
                    # euc dist from all tetra to this pt
                    ndist = np.sqrt(np.sum((tet - up_tet(l)) ** 2, axis=1))
                    if not allow_zero:
                        ndist[np.where(ndist == 0)] = np.inf
                    dist = np.min([dist, ndist], axis=0)
        return dist

    gd = lambda g1, g2: min_dist(rand_genes(g1, M_DRAWS), g2, \
      allow_zero=(set(g1).intersection(g2) == set()))

    #gene_list = ['pr', 'pufM', 'pufL']
    gene_list = all_genes

    gs = gridspec.GridSpec(len(gene_list), len(gene_list))
    gs.update(wspace=0, hspace=0)
    xs = np.linspace(0, 0.1, 200)

    for i, g1 in enumerate(gene_list):
        for j, g2 in enumerate(gene_list):
            ax = plt.subplot(gs[i + j * len(gene_list)])
            ax.text(0.5, 0.95, g1 + '->' + g2, fontsize=2, \
              va='top', ha='center', transform=ax.transAxes)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            dist = gd([g1], [g2])
            ctrl = gd([g1, g2], [g1, g2])
            ax.plot(xs, gaussian_kde(dist)(xs), 'k-')
            ax.plot(xs, gaussian_kde(ctrl)(xs), 'r-')
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
