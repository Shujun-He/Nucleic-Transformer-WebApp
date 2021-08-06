import os

os.environ["ARNIEFILE"] = f"arnie.conf"
import numpy as np
from arnie.bpps import bpps
from decimal import Decimal
import pandas as pd
from tqdm import tqdm
import multiprocessing
from arnie.mea.mea import MEA


def get_predicted_loop_type(sequence, structure, debug=False):
    os.system(f'echo {sequence} > a.dbn')
    os.system('echo \"{}\" >> a.dbn'.format(structure))
    os.system('perl bpRNA/bpRNA.pl a.dbn')
    result = [l.strip('\n') for l in open('a.st')]
    if debug:
        print(sequence)
        print(structure)
        print(result[5])
    return result



def generate_RNA_features(sequence, pkg='rnasoft',Ts=[37,50]):
    bp_matrix_seq=[]
    for T in Ts:
        bp_matrix = bpps(sequence, package=pkg, T=T)
        bp_matrix_seq.append(bp_matrix)
    bp_matrix_seq=np.stack(bp_matrix_seq,0)


    structures=[]
    log_gamma=0
    for j in range(len(bp_matrix_seq)):
        mea_mdl = MEA(bp_matrix_seq[j],gamma=10**log_gamma)
        structures.append(mea_mdl.structure)

    loops=[]
    for j in range(len(structures)):
        # result=get_predicted_loop_type(sequence, structures[j], debug=False)
        # loop=result[5]
        loops.append('X')

    return bp_matrix_seq, structures, loops

#sequence='GGGAAUAAACUAGUAUUCUUCUGGUCCCCACAGACUCAGAGAGAACCCGCCACCAUGAGUAAGGGCGAGGAGCUCUUCACCGGGGUGGUGCCCAUCCUGGUGGAGCUCGACGGGGACGUGAACGGGCACA'
#bp_matrix_seq, structures, loops=generate_RNA_features(sequence)



# plt.imshow(bp_matrix_seq.mean(0))
# plt.show()
# plt.clf()
