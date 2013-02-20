import shutil
import subprocess
import os.path as op
from glob import iglob


def find_classify_orfs(project, scaffolds, markov_model, \
                       mgenemark_loc=None, hmmer_loc=None):
    #TODO: META_MOD should point to actual MetaGeneMark DB
    gene_file = project + '_hmm.txt'
    gff_file = project + '_mgm.gff'
    gfa_file = project + '_mgm.gfa'

    if mgenemark_loc is None:
        mgenemark_loc = ''
    if mgenemark_loc == '':
        from genelinkage.which import which
        for f in which('gmhmmp'):
            mgmloc = op.split(op.realpath(f))[0]
            if op.exists(op.join(mgmloc, 'MetaGeneMark_v1.mod')):
                meta_mod = op.join(mgmloc, 'MetaGeneMark_v1.mod')
                break
        else:
            meta_mod = './MetaGeneMark_v1.mod'
    else:
        meta_mod = op.join(mgenemark_loc, 'MetaGeneMark_v1.mod')
    gmhmmp = op.join(mgenemark_loc, 'gmhmmp')
    aa_from_gff = op.join(mgenemark_loc, 'aa_from_gff.pl')

    subprocess.call([gmhmmp, '-a', '-f', 'G', '-m', meta_mod, '-o',
                     gff_file, scaffolds])
    subprocess.call([aa_from_gff, '<', gff_file, '>', gfa_file])

    if hmmer_loc is None:
        hmmscan = 'hmmscan'
    else:
        hmmscan = op.join(hmmer_loc, 'hmmscan')
    subprocess.call([hmmscan, '--tblout', gene_file, markov_model, gfa_file])


def build_hmmdb():
    GENE_DIR = '/n/pearsonfs1/Roderick/create_COGs/Genes'
    ALL_GENES = GENE_DIR + "all.sto"
    GENE_DB = GENE_DIR + "genes.hmm"
    #TODO: do the blasting on all the individual gene *.FA files
    # to assemble the aligned sequence DBs for each one
    all_gene_file = open(ALL_GENES, 'w')
    for filename in iglob(op.join(GENE_DIR, '*.sto')):
        with open(filename, 'r') as f:
            shutil.copyfileobj(f, ALL_GENES)
    all_gene_file.close()

    subprocess.call(["hmmbuild", GENE_DB, ALL_GENES])
    subprocess.call(["hmmpress", GENE_DB])