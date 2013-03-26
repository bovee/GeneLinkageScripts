import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from genelinkage.calc_dists_and_plot import plot_dist


def gen_figures():
    for i in [1, 4, 16]:
        #print('StatExample' + str(i))
        #gen_fake_clustered('StatExample' + str(i), clust=i)
        print('BioExample' + str(i))
        gen_fake_bio('BioExample' + str(i), clust=i)


def gen_fake(project, dimen=4):
    # randomly distributed genes
    d1 = np.random.random((128, dimen))

    # genes linked to d1
    d2 = np.abs(d1 + np.random.normal(scale=0.01, size=d1.shape))[:64]

    # other randomly distributed genes
    d3 = np.random.random((128, dimen))

    # clustered genes
    d4 = np.abs(np.random.normal(loc=0.5, scale=0.05, size=(64, dimen)))

    # more clustered genes
    d5 = np.abs(np.random.normal(loc=0.2, scale=0.01, size=(32, dimen)))
    d5[:, 0] = np.abs(np.random.normal(loc=0.5, scale=0.01, size=32))

    genes = []
    for d in [d1, d2, d3, d4, d5]:
        genes.append(d / d.sum(axis=1)[:, np.newaxis])
    write_and_plot(project, genes, dimen)


def gen_fake_clustered(project, dimen=4, clust=1):
    genes = []
    for i in range(3):
        if i == 0 or i == 2:
            # randomly distributed genes
            d = np.random.random((128, dimen))
        elif i == 1:
            # genes linked to d1
            d_p = genes[0]
            d = np.abs(d_p + np.random.normal(scale=0.01, \
                                              size=d_p.shape))[32:96]
        elif i == 3:
            # clustered genes
            d = np.abs(np.random.normal(loc=0.5, scale=0.05, size=(64, dimen)))
        elif i == 4:
            # more clustered genes
            d = np.abs(np.random.normal(loc=0.2, scale=0.01, size=(32, dimen)))
            d[:, 0] = np.abs(np.random.normal(loc=0.5, scale=0.01, size=32))

        if i != 1:
            d = d / d.sum(axis=1)[:, np.newaxis]
        if i == 0 or i == 2:
            d[:64] /= clust

        genes.append(d)
    write_and_plot(project, genes, dimen)


def gen_fake_bio(project, dimen=4, clust=1):
    genes = []
    # randomly distributed genes
    rd = np.random.random((256, dimen))
    rd = rd / rd.sum(axis=1)[:, np.newaxis]
    rd[64:192] /= clust
    for i in range(4):
        if i == 0:  # A
            # only present in half of organisms
            d_p = rd[:128]
            d = np.abs(d_p + np.random.normal(scale=0.005, size=d_p.shape))
        elif i == 1:  # B
            # genes linked to A
            d_p = rd[32:96]
            d = np.abs(d_p + np.random.normal(scale=0.005, size=d_p.shape))
        elif i == 2:  # C
            # present in other half of organisms
            d_p = rd[128:]
            d = np.abs(d_p + np.random.normal(scale=0.005, size=d_p.shape))
        elif i == 3:  # D
            # "required" genes
            d = rd
        genes.append(d)
    write_and_plot(project, genes, dimen, has_ctrl=True)


def write_and_plot(project, genes, dimen, has_ctrl=False):
    names = ['A', 'B', 'C', 'D']
    #colors = ['k', 'k', 'grey', 'w']
    colors = ['r', 'b', 'y', 'k']
    #markers = ['s', 'D', 'o', '+']
    markers = ['$\\mathrm{' + n + '}$' for n in names]
    with open(project + '_tetra.txt', 'w') as f:
        f.write('Name,' + ','.join('A' + str(i) for i in range(dimen)))
        for n, g in zip(names, genes):  # loop through each gene
            for t in g:  # go through each "tetranucleotide" for that gene
                f.write('\n' + n + ',' + ','.join(str(i) for i in t))

    if has_ctrl:
        plot_dist(project, bio_ctrl='D')
    else:
        plot_dist(project)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for m, c, g in zip(markers, colors, genes):
        ax.scatter(*g.T[0:3], edgecolor=c, s=30, alpha=1, marker=m)

    # make legend
    #p_styles = [plt.Line2D([0], [0], marker=m, c='w', markerfacecolor=c) \
    #            for c, m in zip(colors, markers)]
    ## readjust first one to make cross not disappear
    #p_styles[-1] = plt.Line2D([0], [0], marker=markers[-1], \
    #                         c='w', markeredgecolor=colors[-1])
    #plt.legend(p_styles, names, numpoints=1)
    plt.show()
