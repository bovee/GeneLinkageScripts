import os.path as op
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.cm import jet as colormap
from matplotlib import mpl
from matplotlib.patches import FancyArrowPatch
import numpy as np


def make_network(folder, projects, cols=None):
    import networkx as nx
    probs = np.array([0.001, 0.0001, 0.00001])
    analz_genes = []
    linked_genes = []

    G = nx.MultiDiGraph()
    for project in projects:
        analz_genes.append(set())
        linked_genes.append(set())
        dist_file = op.join(folder, project + '_tetra_dists.txt')

        for l in open(dist_file, 'r'):
            d = l.split(',')
            prob = float(d[2])  # statistical control
            prob2 = float(d[3])  # biological control
            analz_genes[-1].update(set([d[0], d[1]]))
            if prob < probs[0]:
                linked_genes[-1].update(set([d[0], d[1]]))
                G.add_edge(d[0], d[1], pj=project, pr=prob, pr2=prob2)

    #pos = nx.spring_layout(G, iterations=5)
    pos = nx.circular_layout(G)

    #for prob_type in ['pr', 'pr2']:
    prob_type = 'pr'
    if cols is None:
        cols = int(np.round(np.sqrt(len(projects))))
    gs = gridspec.GridSpec(int(np.ceil(len(projects) / cols)), cols)
    gs.update(wspace=0, hspace=0)

    for p, a_genes, l_genes in zip(projects, analz_genes, linked_genes):
        ax = plt.subplot(gs[projects.index(p)])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)

        for a, c, prob in zip([0.1, 0.3, 1], ['r', 'b', 'g'], probs):
            edges = [(u, v) for u, v, d in G.edges(data=True) \
                    if d['pj'] == p and d[prob_type] < prob]
            #nx.draw_networkx_nodes(G, pos)
            nx.draw_networkx_edges(G, pos, edgelist=edges, \
                                edge_color=c, alpha=a, ax=ax)
        for node in G:
            if node in l_genes:
                font_clr = 'black'
            elif node in a_genes:
                font_clr = 'grey'
            else:
                font_clr = 'white'
            nx.draw_networkx_labels(G, pos, labels={node: node}, ax=ax, \
                                    font_size=8, font_color=font_clr)
        ax.text(0.5, 0.5, p, ha='center', va='center',
                transform=plt.gca().transAxes)
    plt.show()


def text_params_from_ang(a):
    # a goes from 90 down to -270
    offa = 15
    if a <= 90 - offa and a > -90 + offa:
        ha = 'left'
    elif a <= -90 - offa and a > -270 + offa:
        ha = 'right'
    else:
        ha = 'center'
    if a > offa or a <= -180 - offa:
        va = 'bottom'
    elif a <= -offa and a > -120 - offa:
        va = 'top'
    else:
        va = 'center'
    ang = a if a > -90 else a + 180
    return (ha, va, ang)


def make_network_nice(folder, projects, cols=None):
    prob_cut = 0.001

    lnk_genes = []
    anl_genes = []

    for project in projects:
        lnk_genes.append([])
        anl_genes.append(set())
        dist_file = op.join(folder, project + '_tetra_dists.txt')

        for l in open(dist_file, 'r'):
            d = l.split(',')
            prob = float(d[2])  # statistical control
            #prob = float(d[3])  # biological control
            anl_genes[-1].update([d[0], d[1]])
            if prob < prob_cut:
                lnk_genes[-1].append([d[0], d[1], prob])

    if cols is None:
        cols = int(np.round(np.sqrt(len(projects))))

    # all genes compared
    #tot_genes = sorted(set(i for j in anl_genes for i in j))
    # only genes that have links somewhere
    tot_genes = set()
    for l_genes in lnk_genes:
        for g1, g2, p in l_genes:
            tot_genes.update([g1, g2])

    #alternative A:
    #tot_genes = sorted(tot_genes)

    #alternative B: manually order the genes to prevent confusion
    sorted_genes = [i.strip() for i in """
    rpoA,fdhA,soxA-allo,soxA-rhodo,dsrA,
    rpoB,frhB,soxY-allo,soxY-rhodo,dsrB,
    cdsA-allo,cdsA-geo,cdsA-rhodo,cdsA-synn,
    rpoC,soxZ-allo,soxZ-rhodo,ligE
    """.split(',')]
    stot_genes = [g for g in sorted_genes if g in tot_genes]
    tot_genes = stot_genes + sorted(tot_genes.difference(stot_genes))

    adj_f = 2 * np.pi / len(tot_genes)
    g_pos = {g: (np.sin(adj_f * i), np.cos(adj_f * i)) \
             for i, g in enumerate(tot_genes)}
    g_ang = {g: 90 - 360 * i / len(tot_genes) for i, g in \
             enumerate(tot_genes)}

    colormap.set_under(colormap(0))
    colormap.set_over(colormap(colormap.N))
    norm = mpl.colors.LogNorm(vmin=10e-10, vmax=prob_cut)
    p_val_to_c = lambda p: colormap(norm(p))

    fig = plt.figure()
    for proj, grey_labels, l_genes in zip(projects, anl_genes, lnk_genes):
        #ax = plt.subplot(gs[projects.index(proj)])
        ax = fig.add_subplot(int(np.ceil(len(projects) / cols)), cols, \
                             projects.index(proj) + 1)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-1.2, 1.2)
        ax.set_frame_on(False)

        black_labels = set()
        for g1, g2, p in l_genes:
            black_labels.update([g1, g2])

        for g in grey_labels.intersection(tot_genes):
            ha, va, ang = text_params_from_ang(g_ang[g])
            ax.text(*g_pos[g], s=g, size=8, color='grey', va=va, \
                    ha=ha, rotation=ang)
        for g in black_labels:
            ha, va, ang = text_params_from_ang(g_ang[g])
            ax.text(*g_pos[g], s=g, size=8, color='black', va=va, \
                    ha=ha, rotation=ang)

        for g1, g2, p in l_genes:
            clr = p_val_to_c(p)
            arrow = FancyArrowPatch(g_pos[g1], g_pos[g2], color=clr, \
              lw=0.3, arrowstyle='simple,head_length=8,head_width=4')
            ax.add_patch(arrow)

        #TODO: turn off
        ax.text(0.5, 0.5, proj, ha='center', va='center',
                transform=plt.gca().transAxes)
    cax = fig.add_axes([0.33, 0.05, 0.33, 0.010])
    mpl.colorbar.ColorbarBase(cax, cmap=colormap, norm=norm, \
                              orientation='horizontal')
    plt.show()
