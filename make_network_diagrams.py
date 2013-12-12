from genelinkage.network_diagram import make_network_nice

#print('Parameter Tests')
#print('MinGenes 5')
#make_network_nice('../DistProbs/', ['MHL7_1e-5_1000_5', 'MHL7_1e-5_2000_5', 'MHL7_1e-5_5000_5', \
#                                    'MHL7_1e-10_1000_5', 'MHL7_1e-10_2000_5', 'MHL7_1e-10_5000_5', \
#                                    'MHL7_1e-15_1000_5', 'MHL7_1e-15_2000_5', 'MHL7_1e-15_5000_5'])
#print('MinGenes 10')
#make_network_nice('../DistProbs/', ['MHL7_1e-5_1000_10', 'MHL7_1e-5_2000_10', 'MHL7_1e-5_5000_10', \
#                                    'MHL7_1e-10_1000_10', 'MHL7_1e-10_2000_10', 'MHL7_1e-10_5000_10', \
#                                    'MHL7_1e-15_1000_10', 'MHL7_1e-15_2000_10', 'MHL7_1e-15_5000_10'])
#print('MinGenes 25')
#make_network_nice('../DistProbs/', ['MHL7_1e-5_1000_25', 'MHL7_1e-5_2000_25', 'MHL7_1e-5_5000_25', \
#                                    'MHL7_1e-10_1000_25', 'MHL7_1e-10_2000_25', 'MHL7_1e-10_5000_25', \
#                                    'MHL7_1e-15_1000_25', 'MHL7_1e-15_2000_25', 'MHL7_1e-15_5000_25'])


print('Final Figures')
print('e-value + contig length')
#make_network_nice('../FinalProbsWCar/', ['MHL7_1e-5_5000_5', 'MHL7_1e-10_1000_5', \
#                                         'MHL7_1e-10_5000_5', 'MHL7_1e-10_2000_5', \
#                                         'MHL7_1e-15_5000_5', 'MHL7_1e-10_5000_5'])
make_network_nice('../FinalProbsWCar/', ['MHL7_1e-5_5000_5', 'MHL7_1e-10_5000_5', 'MHL7_1e-15_5000_5'], 1, 'fig1.svg')
make_network_nice('../FinalProbsWCar/', ['MHL7_1e-10_1000_5', 'MHL7_1e-10_2000_5', 'MHL7_1e-10_5000_5'], 1, 'fig2.svg')
make_network_nice('../FinalProbsWCar/', ['MHL7_1e-10_5000_5', 'MHL7_1e-10_5000_10', 'MHL7_1e-10_5000_25'], 1, 'fig3.svg')
#make_network_nice('../FinalProbsWCar/', ['MHL7_1e-10_5000_5', 'MHL8_1e-10_5000_5', 'MHLS_1e-10_5000_5'])

print('num genes + samples')
#make_network_nice('../FinalProbsWCar/', ['MHL7_1e-10_5000_5', 'MHL7_1e-10_5000_5', \
#                                         'MHL7_1e-10_5000_10', 'MHL8_1e-10_5000_5', \
#                                         'MHL7_1e-10_5000_25', 'MHLS_1e-10_5000_5'])

