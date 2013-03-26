import os.path as op
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
