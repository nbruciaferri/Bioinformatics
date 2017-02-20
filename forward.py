#Forward algoorithm

# a(transitions) = {(s1,s2):P...}
# e(emissions) = {(s1,e1):P...}
# states = [start, s1, s2, s3....., end]
# F = [[][][][]]

#initialisation part
#----> s = columns
#----> i = rows


#=======================================
def matrix(rows, col):
    M = [[0.0 for i in range(len(rows)+2)]for s in col]
    return M

#=======================================

#initialisation part
def forward (a, e, states, seq):
    F = matrix(seq, states)
    F[0][0] = 1


#iteration part
    for j in range(1, len(seq)+1):
        for i in range(len(states)):
            #print states[i]
            for  w in range(len(states)):
                #print states[w]
                #if states[i] != 'B' and states[w] != 'E':
                F[i][j] = F[i][j] + F[w][j-1] * a.get((states[w], states[i]),0)
                #else:
                    #continue
            #if states[i] != 'B' and states[i] != 'E':
            F[i][j] = F[i][j] * e.get((states[i], seq[j-1]),0)
            #else:
                #continue

#termination part
    for w in range(len(states)):
        F[-1][-1] = F[-1][-1] + F[w][-2] * a.get((states[w], states[i]),0)
    
    return F
    #return F[-1][-1]
#=======================================

#print forward ({('B','Y'):0.2 , ('B','N'):0.8 , ('Y','Y'):0.7 , ('Y','N'):0.2 , ('N','Y'):0.1 , ('N','N'):0.8 , ('Y','E'):0.1 , ('N','E'):0.1 }, {('Y','A'):0.1 , ('Y','G'):0.4 , ('Y','C'):0.4 , ('Y','T'):0.1 , ('N','A'):0.25 , ('N','G'):0.25 , ('N','C'):0.25 , ('N','T'):0.25}, ['B', 'Y', 'N', 'E'], 'ATGCG')
