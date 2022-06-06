import numpy as np


seq1 = str(input('Choose the first sequence to compare: '))
seq2 = str(input('Choose the second sequence to compare: '))

gap_score = int(input('Set te value for the gap score: '))
match_score = int(input('Set te value for the match score: '))
mismatch_score = int(input('Set te value for the mismatch score: '))


def compute_matrix(seq1, seq2):
    '''compute_matrix(seq1,seq2) builds a matrix with an additional row and
    column that represent the increasing gap score'''
    n = len(seq1)
    m = len(seq2)
    gap = gap_score
    match = match_score
    mismatch = mismatch_score

    matrix = np.zeros((n+1,m+1), dtype=int) #builds a matrix with only zeros

    #adding the gap score on the first row and column
    i = 0
    e = 1
    while e < len(matrix[0]): #matrix[0] because I'm considerng the elements of a row
        l = matrix[0]
        l[e] = i + gap
        e += 1
        i += gap
    
    el = 1
    i = 0
    while el < len(matrix): #here I'm considering the columns of the array
        l = matrix
        l[el][0] = i + gap
        el +=1
        i += gap

    return matrix

def NW_matrix(seq1, seq2):
    '''NW_matrix(seq1, seq2) completes the matrix basing on the match, mismatch
    and gap according to the sequences'''
    m = compute_matrix(seq1, seq2)
    gap = gap_score
    match = match_score
    mismatch = mismatch_score
    column = 1
    row = 1    
    while column < len(seq1)+1:
        '''len(seq1)+1 because the counter starts at 1 and not 0'''
        while row < len(seq2)+1:
            v1 = m[column][row-1] + gap
            v2 = m[column-1][row] + gap
            if seq1[column-1] == seq2[row-1]:
                v3 = m[column-1][row-1] + match
                values = [v3,v1, v2]
                m[column][row] = max(values)
                
            if seq1[column-1] != seq2[row-1]:
                v3 = m[column-1][row-1] + mismatch
                values = [v3, v2, v1]
                m[column][row] = max(values)
                
            row += 1
        row = 1
        column += 1
    
    return m
print(NW_matrix(seq1, seq2))
def traceback_matrix(seq1, seq2):
    '''traceback_matrix(seq1, seq2) builds a matrix with only three values:
    1 for diagonal movement, 2 for horizontal movement ( right movement if we
    start from the top of the matrix) and 3 for vertical (down from the top of
    the matrix) movement. Useful for the alignment of the sequences'''
    m = compute_matrix(seq1, seq2)
    gap = gap_score
    match = match_score
    mismatch = mismatch_score
    column = 1
    row = 1
    trbm = np.copy(m) #traceback matrix
    
    #for suppm from the top: 1 = diag, 2 = right, 3 = down

    while column < len(seq1)+1:
        '''len(seq1)+1 because the counter starts at 1 and not 0'''
        while row < len(seq2)+1:
            v1 = m[column][row-1] + gap
            v2 = m[column-1][row] + gap
            if seq1[column-1] == seq2[row-1]:
                v3 = m[column-1][row-1] + match #diagonal
                values = [v3, v1, v2]
                m[column][row] = max(values)
                index = values.index( max(values))
                if index == 0:
                    trbm[column][row] = 1
                if index == 1:
                    trbm[column][row] = 3
                if index == 2:
                    trbm[column][row] = 2
            
            if seq1[column-1] != seq2[row-1]:
                v3 = m[column-1][row-1] + mismatch
                values = [v3, v2, v1]
                m[column][row] = max(values)
                index = values.index( max(values))
                if index == 0:
                    trbm[column][row] = 1
                if index == 1:
                    trbm[column][row] = 3
                if index == 2:
                    trbm[column][row] = 2
            row += 1
        row = 1
        column += 1
    return trbm


def NW_alignment(seq1, seq2):
    '''NW_alignment(seq1, seq2) aligns the two sequences based on the traceback
    matrix'''
    matrix = traceback_matrix(seq1, seq2)
    y = len (matrix)-1
    x = len (matrix[0])-1
    aligned1 = ''
    aligned2 = ''
    #from bottom 1 = diag, 2 = left, 3 = up
    while y >= 0:
        while x >= 0:
            if matrix[y][x] == 1:
                aligned1 += seq1[y-1]
                aligned2 += seq2[x-1]
                y -= 1
                x -= 1
                
            elif matrix[y][x] == 2:
                aligned1 += '-'
                aligned2 += seq2[x-1]
                x -= 1

            elif matrix[y][x] == 3:
                aligned1 += seq1[y-1]
                aligned2 += '-'
                y -= 1
            elif matrix[y][x] <= 0:
                return aligned1[::-1] + '\n' + aligned2[::-1]
                '''as the computed aligned sequences were made backwards,
                they are reversed here and put one below the other in
                order to have a better comparison'''

def NW_traceback():
    '''NW_traceback() builds the traceback over the traceback_matrix()function
    by trasforming the array in list and adding the values D,H and V,
    corresponding to diagonal,horizontal or vertical movement'''
    matrix = traceback_matrix(seq1,seq2)
    matrix = matrix.tolist()
    y = len (matrix)-1 # y as for the y axis ( or column) of the matrix
    x = len (matrix[0])-1 # x for the x axis ( or row) of the matrix
    #from bottom 1 = diagonal, 2 = horizontal, 3 = vertical
    while y >= 0:
        while x >= 0:
            if matrix[y][x] == 1:
                matrix[y][x] = 'D' #D stands for diagonal movement
                y -= 1
                x -= 1
                
            elif matrix[y][x] == 2:
                matrix[y][x] = 'H' #H stands for horizontal movement
                x -= 1

            elif matrix[y][x] == 3:
                matrix[y][x] = 'V' #V stands for vertical
                y -= 1
                
            elif matrix[y][x] <= 0: #it stops because it has found elemements 
                return matrix       #from first row or first column

print(NW_alignment(seq1, seq2))
print(NW_traceback())

