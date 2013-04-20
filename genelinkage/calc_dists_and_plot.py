# need to do this first to make plots in file
import matplotlib
matplotlib.use('Agg')

import os.path as op
from collections import Counter
from functools import partial
from multiprocessing import Pool

import numpy as np
import scipy.spatial.distance as ssd
from scipy.stats import gaussian_kde  # , mannwhitneyu
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

MULTIPROCESSING = False

try:
    from rpy2.robjects import r, FloatVector
    ranksum = lambda a, b: r['wilcox.test'](FloatVector(a), FloatVector(b), \
                           alternative='less').rx('p.value')[0][0]
except:
    ranksum = None
    print('RPy not loaded. No stat calc done.')


def min_dst(tet1, tet2, allow_zero=True):
    """
    Given two sets of arrays, find the minimum distance from
    each member of tet1 to its nearest neighbor in tet2
    """
    dists = ssd.cdist(tet1, tet2)
    if not allow_zero:
        dists[dists == 0] = np.inf
    return dists.min(axis=1)

    #dists = np.empty(tet1.shape[0])
    #for i, t1 in enumerate(tet1):
    #    min_dist = np.sum((tet2 - t1) ** 2, axis=1)
    #    if not allow_zero:
    #        dists[i] = np.min(min_dist[min_dist != 0])
    #    else:
    #        dists[i] = np.min(min_dist)
    #return np.sqrt(dists)


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


def samp_dist(gene2, gene1, gene_data=None, ctrl=False, monte_draws=None):
    """
    Given two gene names, find the minimum distances from 1 to 2.
    """
    gene_tetra, gene_ct, gene_ids, gene_names = gene_data
    if ctrl:
        gene_locs = np.logical_or(gene_ids == gene_names[gene1], \
                                  gene_ids == gene_names[gene2])
        dist = []
        for _ in range(monte_draws):
            rand_picks = np.random.permutation(sum(gene_locs))
            tet1 = gene_tetra[gene_locs][rand_picks[:gene_ct[gene1]]]
            tet2 = gene_tetra[gene_locs][rand_picks[-gene_ct[gene2] - 1:]]
            dist += list(min_dst(tet1, tet2, allow_zero=False))
        #tet1 = gene_tetra[gene_ids == gene_names[gene1]]
        #dist = []
        #for _ in range(monte_draws):
        #    tet2 = rnd_genes([], gene_ct[gene1], gene_data)
        #    dist += list(min_dst(tet1, tet2, allow_zero=(gene1 != gene2)))
    else:
        tet1 = gene_tetra[gene_ids == gene_names[gene1]]
        tet2 = gene_tetra[gene_ids == gene_names[gene2]]
        dist = list(min_dst(tet1, tet2, allow_zero=(gene1 != gene2)))
    return dist


def plot_dist(project, gene_list=None, filt_length=None, monte_draws=100, \
              plot=True, save=False, bio_ctrl='rpoA', save_file=None):
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
        gene_list = sorted(list(gene_ct.keys()))

    if filt_length is None:
        gene_list = [g for g in gene_list if g in gene_ct]
    else:
        gene_list = [g for g in gene_list if gene_ct[g] >= filt_length \
                     and g in gene_ct]

    if MULTIPROCESSING:
        po = Pool()

    if gene_ct.get(bio_ctrl, 0) > 5:
        #if bio_ctrl in gene_list:
        #    del gene_list[gene_list.index(bio_ctrl)]
        ctrl_rpo = samp_dist(bio_ctrl, bio_ctrl, gene_data, \
                                monte_draws=monte_draws)
    else:
        ctrl_rpo = None

    if save:
        fname = op.splitext(op.realpath(tetra_file))[0]
        if save_file is not None:
            fname += '_' + save_file
        else:
            fname += '_dists.txt'
        dist_file = open(fname, 'w')
    if plot:
        gs = gridspec.GridSpec(len(gene_list), len(gene_list))
        gs.update(wspace=0, hspace=0)
        xs = np.linspace(0, 0.1, 200)
        if hasattr(gaussian_kde, 'set_bandwidth'):
            kdeargs = [0.1]
        else:
            kdeargs = []

    for i, g1 in enumerate(gene_list):
        print(int(100 * i / len(gene_list)), g1)
        dist_f = partial(samp_dist, gene1=g1, gene_data=gene_data, \
                         monte_draws=monte_draws)
        ctrl_f = partial(samp_dist, gene1=g1, gene_data=gene_data, \
                         ctrl=True, monte_draws=monte_draws)
        if MULTIPROCESSING:
            dists = po.map(dist_f, gene_list)
            ctrls = po.map(ctrl_f, gene_list)
        else:
            dists = list(map(dist_f, gene_list))
            ctrls = list(map(ctrl_f, gene_list))
        for j, g2 in enumerate(gene_list):
            if ranksum is not None:
                #mwu = mannwhitneyu(ctrls[j], dists[j])[1]
                #mwu2 = mannwhitneyu(ctrl_rpo, dists[j])[1]
                mwu = ranksum(dists[j], ctrls[j])
                if ctrl_rpo is not None:
                    mwu2 = ranksum(dists[j], ctrl_rpo)
                else:
                    mwu2 = 0
            else:
                mwu, mwu2 = None, None

            if save and mwu is not None:
                dist_file.write(g1 + ',' + g2 + ',' + \
                                str(mwu) + ',' + str(mwu2) + '\n')

            if plot:
                ax = plt.subplot(gs[i + j * len(gene_list)])
                if mwu is not None and g1 != g2:
                    txt = g1 + '$\\rightarrow$' + g2 + \
                        '\np$_{s}$=' + '{:.2%}'.format(mwu) + \
                        ' (p$_{b}$=' + '{:.2%})'.format(mwu2)
                elif mwu2 is not None:
                    txt = g1 + '$\\rightarrow$' + g2 + \
                        '\np$_{b}$=' + '{:.2%}'.format(mwu2)
                else:
                    txt = g1 + '$\\rightarrow$' + g2
                ax.text(0.5, 0.95, txt, fontsize=6, \
                        va='top', ha='center', transform=ax.transAxes)
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                if sum(ctrls[j]) != 0:
                    ys = gaussian_kde(ctrls[j], *kdeargs)(xs)
                    ax.plot(xs, ys, '-', c='0.3')
                if sum(dists[j]) != 0:
                    ys = gaussian_kde(dists[j], *kdeargs)(xs)
                    ax.plot(xs, ys, 'k-')
                if ctrl_rpo is not None:
                    ys = gaussian_kde(ctrl_rpo, *kdeargs)(xs)
                    ax.plot(xs, ys, '--', c='0.3')
                    ax.set_ybound(0, 1.2 * np.max(ys))
    if save:
        dist_file.close()
    if plot:
        plt.gcf().set_size_inches(24, 24)
        fig_file = op.splitext(op.realpath(tetra_file))[0] + '.png'
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')

    del gene_data, gene_tetra, gene_ct, gene_ids, gene_names

if __name__ == '__main__':
    # run sensitivity tests on the monte carlo
    proj = '/n/home04/bovee/PLFolder/GeneLinkage/MHLParams/MahoneyLake7M_10e-10_2000_10'
    for mtc in [100, 500, 1000, 5000]:
        for i in range(5):
            plot_dist(proj, filt_length=5, monte_draws=mtc, plot=False, \
                      save=True, save_file='test_' + str(i) + '_' + str(mtc) \
                      + '.txt')
