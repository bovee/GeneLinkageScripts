#!/bin/bash
module load bio/ncbi-blast-2.2.25+
BLAST_DB=/n/pearsonfs1/SequenceData/NRBlastDB/nr
DATA_DIR=/n/pearsonfs1/Roderick/create_COGs/Genes
GENE=$1

psiblast -db ${BLAST_DB} -query ${DATA_DIR}/${GENE}.fa -outfmt "6 sseqid" -out ${DATA_DIR}/${GENE}.blastout -evalue 1e-10 -num_alignments 100000 -num_descriptions 100000 -max_target_seqs 100000
blastdbcmd -db ${BLAST_DB} -entry_batch ${DATA_DIR}/${GENE}.blastout -out ${DATA_DIR}/${GENE}.blast.fa

module load bio/muscle-3.8.31
muscle -in ${DATA_DIR}/${GENE}.blast.fa -out ${DATA_DIR}/${GENE}.aln.fa -seqtype protein -maxiters 2

module load hpc/biopython-1.60_python-2.7.3
python -c "import Bio.AlignIO; Bio.AlignIO.convert('${DATA_DIR}/${GENE}.aln.fa','fasta','${DATA_DIR}/${GENE}.sto','stockholm')"
sed -i '1 a \#=GF ID $GENE' ${DATA_DIR}/${GENE}.sto
