import itertools


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


def generate_tetra(project, scaffolds, min_len=2000):
    """
    gene_file: HMMER predicted genes.

    """
    #TODO: screed is only for python 2, so can't use this on py3
    import screed

    gene_file = project + '_hmm.txt'
    tetra_file = project + '_tetra.txt'
    gff_file = project + '_mgm.gff'

    # generate a mapping to merge reverse complement tetranucleotides
    seq_map = {'': 4 * 'N'}
    for s in (''.join(i) for i in itertools.product(*(4 * ['ATGC']))):
        if rc(s) not in seq_map or s not in seq_map:
            seq_map[s] = s
            seq_map[rc(s)] = s

    # write out the headers for the gene/tetranucleotide file
    fout = open(tetra_file, 'w')
    srted_vals = list(set(seq_map.values()))
    srted_vals.sort()
    fout.write(','.join(['Gene'] + srted_vals) + '\n')

    f = open(gff_file, "r")
    for i in range(7):
        f.readline()

    # make a dictionary to map genes back to their contigs
    gene2contig = {}
    for ln in f:
        flds = ln.strip().split('\t')
        gene2contig[flds[-1].replace(' ', '_')] = flds[0].split(' ')[0]
    f.close()

    hmm = open(gene_file, 'r')
    contigs = screed.read_fasta_sequences(scaffolds)

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
            cc = contigs[gene2contig[gff_name]]
            if len(cc) < min_len:
                continue
            frq = dict([(s, 0) for s in seq_map.values()])
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
