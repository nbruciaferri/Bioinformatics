#####################################################################
#| This code was created by Niccolo' Bruciaferri.                  |#
#| It implements the HMM Viterbi algorithm for alignment of        |#
#| biological sequences in python.                                 |#
#|                                                                 |#
#|                                                                 |#
#|                                                                 |#
#|                                                                 |#
#| Master degree in Bioinformatics                                 |#
#| Bologna, 2017-02-08                                             |#
#####################################################################
#-------------------------------------------------------------------|
#Definition of the variables

a = {('B','Y'):0.2 , ('B','N'):0.8 , ('Y','Y'):0.7 , ('Y','N'):0.2 , ('N','Y'):0.1 , ('N','N'):0.8 , ('Y','E'):0.1 , ('N','E'):0.1 }
#a is a dictionary of transitions

e = {('Y','A'):0.1 , ('Y','G'):0.4 , ('Y','C'):0.4 , ('Y','T'):0.1 , ('N','A'):0.25 , ('N','G'):0.25 , ('N','C'):0.25 , ('N','T'):0.25}
#e is a dictionary of emissions

states = ['B', 'Y', 'N', 'E']
#list of states 

seq = 'ATGCG'
#target sequence
#-------------------------------------------------------------------|
#====================================================================
def viterbi(a,e,states,seq):

    '''The principal approach used in this code in order to implement
   the Viterbi algorithm, is to store the maximum in lists, in 
   order to keep the traceback of where the maximum is taken'''

#Initialisation
    V = [[0.0 for i in range(len(seq)+2)]for s in states]
    V[0][0] = 1.0

    B = [[0.0 for i in range(len(seq)+2)]for s in states]
    #B is the backtrace matrix
    path = ''

    #for i in range(1,len(states)):
        #V[i][0] = a.get((states[0], states[i]),0.0)
        #Viterbi path for j = 0

#Recurrence
    for j in range(1,len(seq)+1):
    #Viterbi path for j > 0
        for i in range(len(states)):
            b = []
            #b is used as accumulator
            n = 0
            for w in range(len(states)):
                #V[i][j] = max([V[w][j-1]*a.get((states[w], states[i]),0.0)*e.get((states[i], seq[j-1]),0.0)])
                n = V[w][j-1] * a.get((states[w], states[i]),0.0)
                b.append(n)
            tmp = max(b)
            for w in range(len(states)):
                if tmp == b[w]:
                    B[i][j] = states[w]
            V[i][j] = tmp * e.get((states[i], seq[j-1]),0.0)

#Termination
    T = []
    for w in range(len(states)):
        #V[-1][-1] = max([V[w][-2]*a.get((states[w], states[-1]),0.0)])
        T.append(V[w][-2]  * a.get((states[w], states[-1]),0.0))
    V[-1][-1] = max(T)
    for k in range(len(states)):
        if V[-1][-1] == T[k]:
            B[-1][-1] = states[k]
#backtrace
    tmp = 0
    for j in range(1,len(seq)+1):
        N = []
        for i in range(len(states)):
            N.append(V[i][j])
        tmp = max(N)
        for i in range(len(states)):
            if N[i] == tmp:
                path = path + B[i][j]

    #for i in range(1,len(seq)+1):
        #l = []
        #for j in range(1,len(states)-1):
            #l.append(V[j][i])
            #print l
        #for j in range(1,len(states)-1):
            #if V[j][i]  == max(l):
                #path = path + states[j]
    #print B
    return path + 'E'
    #return V
#====================================================================

print viterbi(a,e,states,seq)
