import itertools

def perm(s):

    '''perm(s) creates the possible permutations between elements in s (in this case
    the k-mers).'''
    
    p = list(itertools.permutations(s))
    return p #a list containing tuples


def overlap(e): #e is an element from perm(s), a tuple

    '''overlap(e) merge or overlap the k-mers into a single string'''
    
    e = list(e)
    while len(e) != 1:
        i = 0
        a = len(e[0])
        b = len(e[1])
        counter = min(a,b)
            
        while counter > 0:
            if e[0][i:] == e[1][:counter]:
                e[0] = e[0][0:counter+2] + e[1][counter:]
                print(e[0])
                e.remove(e[1])
                counter = -1
            else:
                i +=1
                counter -= 1
                if counter == 0:
                    e[0] = e[0] + e[1]
                    e.remove(e[1])
                    counter = -1
        
    return e

def bruteforce(s):

    '''bruteforce(s) does the overlap for each permutation and then selects the
    shortest superstrings that it gets'''
    
    p = perm(s)
    all_sol = [overlap(e) for e in p]
    sol = []
    for sublist in all_sol:
        for item in sublist:
            sol.append(item)
    
    length = len(min(sol, key = len))
    shortest = []
    for e in sol:
        if len(e) == length:
            shortest.append(e)
    return shortest

if name == '__main__':
    if len(sys.argv) != 2:
        print("Ussage: %s list of k-mers s " % sys.argv[0])
        sys.exit()
    s = sys.argv[1]
    bruteforce(s)
        
                
                
