# need to do this first to make plots in file
import matplotlib
matplotlib.use('Agg')

import os.path as op
from collections import Counter
from functools import partial
from multiprocessing import Pool

import numpy as np
from scipy.stats import gaussian_kde  # , mannwhitneyu
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

try:
    from rpy2.robjects import r, FloatVector
    ranksum = lambda a, b: r['wilcox.test'](FloatVector(a), FloatVector(b), \
              alternative='greater').rx('p.value')[0][0]
except:
    ranksum = None


def min_dst(tet1, tet2, allow_zero=True):
    """
    Given two sets of arrays, find the minimum distance from
    each member of tet1 to its nearest neighbor in tet2
    """
    dists = np.empty(tet1.shape[0])
    for i, t1 in enumerate(tet1):
        min_dist = np.sum((tet2 - t1) ** 2, axis=1)
        if not allow_zero:
            dists[i] = np.min(min_dist[min_dist != 0])
        else:
            dists[i] = np.min(min_dist)
    return dists


def rnd_genes(genes=[], n=1, gene_data=None):
    """
    Pick n random HMMER-identified "genes" with replacement
    from the gene-tetranucleotide file.

    If genes == [] pick from all genes.
    """
    if gene_data is None:
        return np.array([])
    gene_tetra, gene_ct, gene_ids, gene_names = gene_data
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


def samp_dist(gene2, gene1, gene_data=None, monte_draws=None):
    """
    Given two gene names, find the minimum distances from 1 to 2.
    """
    gene_tetra, gene_ct, gene_ids, gene_names = gene_data
    tet1 = gene_tetra[gene_ids == gene_names[gene1]]
    if gene2 is None:
        dist = []
        for _ in range(monte_draws):
            tet2 = rnd_genes([], gene_ct[gene1], gene_data)
            dist += list(min_dst(tet1, tet2, allow_zero=(gene1 != gene2)))
    else:
        tet2 = gene_tetra[gene_ids == gene_names[gene2]]
        dist = list(min_dst(tet1, tet2, allow_zero=(gene1 != gene2)))
    return dist


def plot_dist(project, gene_list=None, filt_length=None, monte_draws=10, \
              plot=True, save=False):
    """
    gene_list: list of genes to compare. If None, compare all genes
        listed in the tetra file
    filt_length: genes who don't have at least this number of
        ORFs are not compared
    monte_draws: number of Monte Carlo draws to do for controls
    """
    tetra_file = project + '_tetra.txt'

    # convenience functions to return a gene name or
    # its tetranucleotide array
    up_name = lambda l: l.split(',')[0].strip()
    up_tet = lambda l: np.array([float(i) for i in l.split(',')[1:]])

    print('Buffering')
    # first, get a count of how many of each gene there are
    with open(tetra_file, 'r') as f:
        # how many dimensions are in each tetra score?
        tetlen = len(f.readline().split(',')) - 1
        gene_names, gene_ids = {}, []
        # get a count of how many of each gene there are
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
    gene_data = (gene_tetra, gene_ct, gene_ids, gene_names)
    print('Done buffering')

    if gene_list is None:
        gene_list = list(gene_ct.keys())
    if filt_length is None:
        gene_list = [g for g in gene_list if g in gene_ct]
    else:
        gene_list = [g for g in gene_list if gene_ct[g] >= filt_length \
                     and g in gene_ct]

    po = Pool()

    if save:
        fname = op.splitext(op.realpath(tetra_file))[0] + '_dists.txt'
        dist_file = open(fname, 'w')
        #mk_str = lambda d: ','.join([str(i) for i in d])
    if plot:
        gs = gridspec.GridSpec(len(gene_list), len(gene_list))
        gs.update(wspace=0, hspace=0)
        xs = np.linspace(0, 0.005, 200)

    ctrl_rpo = samp_dist('rpoA', 'rpoA', gene_data, \
                         monte_draws=monte_draws)

    #if save:
    #    dist_file.write('rpo_ctrl,' + mk_str(ctrl_rpo))

    for i, g1 in enumerate(gene_list):
        print(int(100 * i / len(gene_list)), g1)
        dist_f = partial(samp_dist, gene1=g1, gene_data=gene_data, \
                         monte_draws=monte_draws)
        dists = po.map(dist_f, gene_list)
        ctrls = po.map(dist_f, [None] * len(gene_list))
        #dists = list(map(dist_f, gene_list))
        #ctrls = list(map(dist_f, [None] * len(gene_list)))
        for j, g2 in enumerate(gene_list):
            if ranksum is not None:
                #mwu = mannwhitneyu(ctrls[j], dists[j])[1]
                #mwu2 = mannwhitneyu(ctrl_rpo, dists[j])[1]
                mwu = ranksum(dists[j], ctrls[j])
                mwu2 = ranksum(dists[j], ctrl_rpo)
            else:
                mwu, mwu2 = None, None

            if save and mwu is not None:
                dist_file.write(g1 + ',' + g2 + ',' + \
                                str(mwu) + ',' + str(mwu2))
                #dist_file.write(g1 + '->' + g2 + ',' + mk_str(dists[j]))
                #dist_file.write(g1 + '->' + g2 + '_ctrl,' + mk_str(ctrls[j]))

            if plot:
                ax = plt.subplot(gs[i + j * len(gene_list)])
                if mwu is not None:
                    txt = g1 + '->' + g2 + \
                            '\np={:.2e} ({:.2e})'.format(mwu, mwu2)
                else:
                    txt = g1 + '->' + g2
                ax.text(0.5, 0.95, txt, fontsize=2, \
                va='top', ha='center', transform=ax.transAxes)
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                ax.plot(xs, gaussian_kde(ctrl_rpo, 0.1)(xs), 'b-')
                if sum(dists[j]) != 0:
                    ax.plot(xs, gaussian_kde(ctrls[j], 0.1)(xs), 'r-')
                if sum(dists[j]) != 0:
                    ax.plot(xs, gaussian_kde(dists[j], 0.1)(xs), 'k-')
    if save:
        dist_file.close()
    if plot:
        plt.gcf().set_size_inches(24, 24)
        fig_file = op.splitext(op.realpath(tetra_file))[0] + '.png'
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')
