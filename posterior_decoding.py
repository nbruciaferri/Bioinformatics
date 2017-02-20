#####################################################################
#|  This code was created by Niccolo' Bruciaferri.                 |#
#|  It implements the posterior decoding algorithm in python       |#
#|  language.                                                      |#
#|                                                                 |#
#|  Master Degree in Bioinformatics,                               |#                       
#|  Bologna, 2017-02-11                                            |#
#####################################################################

#Importing forward and backword matrices
from forward import *
from backward import *
#-------------------------------------------------------------------|
#Definition of the variables
F = forward(a,e,states,seq)
B = backward(a,e,states,seq)

states = ['B', 'Y', 'N', 'E']

seq = 'ATGCG'
#-------------------------------------------------------------------|
#====================================================================
def posterior_decoding(seq, states):
    PD = matrix(seq, states)
    BEST_PATH_PD = ''
    #print F
    #print B
    for j in range(1,len(seq)+1):
        #print seq[j]
        for i in range(len(states)):
            PD[i][j] =(F[i][j]*B[i][j])/F[-1][-1]
    #print PD
    for j in range(1,len(seq)+1):
        (PD[i][j], state) = max((PD[i][j], states[i])for i in range(len(states)))
        #print state
        BEST_PATH_PD = BEST_PATH_PD + state

    return BEST_PATH_PD
#====================================================================
print posterior_decoding(seq, states)
