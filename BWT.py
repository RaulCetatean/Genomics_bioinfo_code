import numpy as np
import sys

def BWT_l(T):

    '''BWT takes a string T, adds the $ sign and then computes all the cyclic rotations
    of that string as a list of lists'''

    T = T + '$'
    matrix = []

    for e in range(len(T)):
        ns = T[e+1:] + T[:e+1]
        matrix = matrix + [ns]
    return matrix

def BWT_array():

    '''BWT_array takes the output of BWT, sorts it and transforms it into an array for
    better visual representaion'''

    matrix =sorted(BWT_l(T))
    matrix = [list(e) for e in matrix]
    array = np.asarray(matrix)
    return array
    
def L_col():

    '''L_col considers only the last column of the sorted BWT'''

    array = BWT_array()
    L =[]
    v = len(array[0])

    for e in array:
        L = L + [e[v-1]]
    return L

def F_col():

    '''F_col considers only the first column of the sorted BWT'''

    F = sorted(L_col())
    return F
    
def rev_BWT():

    '''rev_BWT reconstructs the original string s from the last column of the sorted BWT'''

    L = L_col()
    F = F_col()
    A = []
    rev_s = ''

    for e in L:
        A += [[e, F.index(e)]]
        F[F.index(e)] = 0
        '''The value is set to zero in order to get the indexes of the other
        characters that are the same but at different index. This helps to go
        from the character present in the last column to the same character in
        the first column (characters with the same rank)'''
        
    i = 0
    while i >= 0:
        rev_s += A[i][0]
        i = A[i][1]
        if rev_s[-1] == '$':
            i = -1  #value set to -1 in order to exit the loop
    return rev_s

def query(P):

	'''query controls if the query string P is present in the bwt considering only the first 
	and last column of the transform. Return value is a boolean and the list of indexes (if more
	than one) where P is found (if found)'''

    L = L_col()
    F = F_col()
    f = F[:]
    A = []
    for e in L:
        A += [[e, f.index(e)]]
        f[f.index(e)] = 0
    '''A is the list containing the characters from L and their respective
    index from the F column'''

    new_ind = []
    ind_f = []
    i = len(P)-1
    while i >= 0:
        if len(new_ind) != 0 :
            ind_f = new_ind[:]
        else:
            ind_f = [n for n,val in enumerate(F) if val == P[i]]

        j = 0
        while j < len(ind_f):
            if L[ind_f[j]] != P[i-1]:
                ind_f.remove(ind_f[j])
            elif len(ind_f) == 0:
                return False
            else:
                break
        '''this while is a control to see if the next character if found in the L column.
        If is not found, the list becomes empty and returns False'''

        i -= 1
        ind_l = [n for n,val in enumerate(L) if val == P[i]]
            
        
        ind = [e for e in ind_f if e in ind_l]
        for e in ind:
            new_ind += A[e]
        new_ind = new_ind[1::2]
        '''since in new_ind I have always first a character and then a number, I slice
        new_ind in order to consider only the numbers(that are the indexes'''
        #now I have the indexes that should go to F column
        k = 0
        while j < len(new_ind):
            if F[new_ind[k]] != P[i-1]:
                new_ind.remove(new_ind[k])
            elif len(new_ind) == 0:
                return False
            else:
                break
        '''same control as previous one but controls if it can go from L column to F column'''
     
        i -= 1
        
    return True, ind_f
    '''result is a tuple with a boolean and index or indexes of F column where
    the query P is found'''
        
def offset(P):

    '''offset gives you the offset for the query string P everytime the string is found in the bwt'''

    ind = query(P) 
    matrix = sorted(BWT_l(s))
    print(matrix)
    r_array = [[len(e)- e.index('$')-1] for e in matrix]
    result = [r_array[e] for e in ind[1]] #a list with the offsets ( in case P is present more than once)
    return result

def BWT(T,P):
	bw_transform = BWT_array()
	query = query(P)
	offset = offset(P)
	return bw_transform, query, offset

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("Usage: %s string T query P" % sys.argv[0] )
		sys.exit()
	T = sys.argv[1]
	S = sys.argv[2]
	BWT(T,P)



