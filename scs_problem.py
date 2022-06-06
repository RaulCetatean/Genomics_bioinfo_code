import sys
import itertools
import random

def combination(s):
    
    '''combinations(s) creates all the combinations of two k-mers from s'''
    
    comb = list(itertools.combinations(s,2))
    return comb #a list of tuples

def overlap():

    '''the function overlap(s) computes the overlap between the couples formed
    with combinations(s) and how many characters overlap.'''
    
    comb = combination(s)
    storage = {}
        
    for e in comb:
        i = 0
        a = len(e[0])
        b = len(e[1])
        counter = min(a,b) #max possible value for the overlap
                         
        while counter > 0:
            if e[0][i:] == e[1][:counter]:
                storage[e[0], e[1], e[0][0:counter+2] + e[1][counter:] ] = counter
                counter = -1 
            else:
                i += 1
                counter -= 1
                if counter == 0:
                    storage[e[0], e[1], e[0]+e[1]] = counter
    '''since the order of the k-mers matters, this second part has the same operations
    but considering first the k-mer previously considered second'''

    for e in comb:
        i = 0
        a = len(e[0])
        b = len(e[1])
        counter = min(a,b)
                
        while counter > 0:
            if e[1][i:] == e[0][:counter]:
                storage[e[1], e[0], e[1][0:counter+2] + e[0][counter:] ] = counter
                counter = -1
            else:                
                i += 1
                counter -= 1
                if counter == 0:
                    storage[e[1], e[0], e[1]+e[0]] = counter
                    counter = -1
                                        
    return storage
print(overlap())
def scs(s):

    '''scs(s,t) iterates through s until it remains only one element, that is the
    product of the merging of the k-mers.'''
    
    while len(s) != 1:
        storage = overlap()
        val = list(storage.values())
        counter = max(val)
        l = []
        nl = []
        for e in storage:
            if storage[e] == counter:
                l += e
                
        nl = [tuple(l[i:i+3]) for i in range(0, len(l), 3)]
                
        if len(nl) == 1:
            s.remove(nl[0][0])
            s.remove(nl[0][1])
            s.insert(0, nl[0][2])
        else:
            n = random.randint(0, len(nl)-1)
            el = nl[n] 
            s.remove(el[0])
            s.remove(el[1])
            s.insert(0, el[2])
    return s[0]

if name == "__main__":
    if len(sys.argv) != 2:
        print("Ussage: %s list of k-mers s " % sys.argv[0])
        sys.exit()
    s = sys.argv[1]
    scs(s)

    




