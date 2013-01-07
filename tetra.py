#/usr/bin/env python
import os
import shutil
import math
import itertools
import subprocess
from math import sqrt
import screed

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

build_hmmdb = False
build_tetra = False
build_dist = True

if build_hmmdb:
    all_gene_file = open(ALL_GENES, 'w')
    for filename in iglob(os.path.join(GENE_DIR, '*.sto')):
        with open(filename, 'r') as f:
            shutil.copyfileobj(f, ALL_GENES)
    all_gene_file.close()
     
    subprocess.call(["hmmbuild", GENE_DB, ALL_GENES])
    subprocess.call(["hmmpress", GENE_DB])
    subprocess.call(["hmmscan", "--tblout", GENE_FILE, GENE_DB, GFA_FILE]])

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


    seq_map = {}
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
        gene2contig[flds[-1].replace(' ','_')] = flds[0].split(' ')[0]
    f.close()

    hmm = open(GENE_FILE, 'r')
    contigs = screed.read_fasta_sequences(FASTA_FILE)

    for i in range(3):
        hmm.readline()

    for ln in hmm:
        gene_name = ln[0:21].strip()

        frq = dict([(s, 0) for s in seq_map.values()])
        cc = contigs[gene2contig[ln[32:53].strip()]]
        for ss in slid_win(str(cc.sequence).upper(), 4):
            frq[seq_map.get(ss, 'NNNN')] += 1

        sum_frq = float(sum(frq.values()))
        if sum_frq == 0:
            sum_frq = 1
        fout.write(','.join([gene_name] + [str(frq[i] / sum_frq) for i in srted_vals]))
        fout.write('\n')
        fout.flush()
    hmm.close()
    fout.close()

# start computing gene distances
if build_dist:
    gene_tetra = open(TETRA_FILE, "r")
    gene_tetra.readline()

    gene_dist = open(DIST_FILE, "w")
    gene_dist.write("FromGene,ToGene,Dist\n")

    for ln in gene_tetra:
        l = ln.split(',')
        gene_tetra_2 = open(TETRA_FILE, "r")
        gene_tetra_2.readline()
        dist_dict = {}
        for ln2 in gene_tetra_2:
            l2 = ln2.split(',')
            dist = math.sqrt(sum((float(i) - float(j))**2 for i,j in zip(l[1:], l2[1:])))
            if not (l[0] == l2[0] and dist == 0):
                dist_dict[l2[0]] = min(dist, dist_dict.get(l2[0], 100))
            for gene in dist_dict:
                gene_dist.write(l[0] + ',' + gene + ',' + str(dist_dict[gene])+'\n')
        gene_dist.flush()
    gene_tetra.close()
    gene_dist.close()

# library(ggplot2)
# d <- read.csv('gene_dist.txt')
# genes <- c('dsrA','dsrB','dsrC','dsrL','dsrN')
# ggplot(subset(d, FromGene %in% c(genes) & ToGene %in% c(genes)), aes(x=Dist)) + facet_grid(FromGene ~ ToGene, labeller=label_both) + geom_density() + xlim(0,1)
