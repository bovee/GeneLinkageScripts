module load bio/hmmer-3.0
DATA_DIR=/n/pearsonfs1/Roderick/create_COGs/Genes

#TODO: make sure all *.sto files end with //
cat ${DATA_DIR}/*.sto > ${DATA_DIR}/all.sto
hmmbuild ${DATA_DIR}/genes.hmm ${DATA_DIR}/all.sto
hmmpress ${DATA_DIR}/genes.hmm
hmmscan --tblout ./hmm_results.txt ${DATA_DIR}/genes.hmm ~/MHL7/Genes/genes.faa > /dev/null
