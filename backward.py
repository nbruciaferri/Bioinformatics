#Backward algorithim

# a(transitions) = {(s1,s2):P...}
# e(emissions) = {(s1,e1):P...}
# states = [start, s1, s2, s3....., end]
# F = [[][][][]]

#----> s = columns
#----> i = rows
#====================================================================
#------> definition of the variables

a ={('B','Y'):0.2 , ('B','N'):0.8 , ('Y','Y'):0.7 , ('Y','N'):0.2 , ('N','Y'):0.1 , ('N','N'):0.8 , ('Y','E'):0.1 , ('N','E'):0.1 }

e = {('Y','A'):0.1 , ('Y','G'):0.4 , ('Y','C'):0.4 , ('Y','T'):0.1 , ('N','A'):0.25 , ('N','G'):0.25 , ('N','C'):0.25 , ('N','T'):0.25}

states = ['B', 'Y', 'N', 'E']

seq = 'ATGCG'

#----->
#====================================================================
def backward(a,e,states,seq):
    B = [[0.0 for i in range(len(seq)+2)]for j in states]
    B[-1][-1] = 1.0
    
    for i in range(len(states)):
        B[i][-2] = a.get((states[i], states[-1]),0)   
    
    for j in range(len(seq)-1,0,-1):
        #print seq[j-1]
        for i in range(len(states)-2,0,-1):
            for w in range(len(states)-2,0,-1):
                B[i][j] = B[i][j] + B[w][j+1] * a.get((states[i], states[w]),0) * e.get((states[w], seq[j]),0)
 
    for w in range(len(states)):
        B[0][0] = B[0][0] + B[w][1] * a.get((states[0], states[w]),0) * e.get((states[w], seq[0]),0)
        #print (states[w], states[-1])
        #print a.get((states[w], states[-1]),0)
        #print (states[w], states[1])
        #print a.get((states[-1], states[w]),0)
    #return B[0][0]
    return B
#====================================================================

#print backward(a,e,states,seq)
