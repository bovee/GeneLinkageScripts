import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def gen_fake(project, dimen=4):
    genes = []
    colors = ['r', 'g', 'b', 'y', 'k']

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

    for d in [d1, d2, d3, d4, d5]:
        genes.append(d / d.sum(axis=1)[:, np.newaxis])

    with open(project + '_tetra.txt', 'w') as f:
        f.write('Name,' + ','.join('A' + str(i) for i in range(dimen)))
        for c, g in zip(colors, genes):  # loop through each gene
            for t in g:  # go through each "tetranucleotide" for that gene
                f.write('\n' + c + ',' + ','.join(str(i) for i in t))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for c, g in zip(colors, genes):
        ax.scatter(*g.T[0:3], c=c)
    plt.show()
