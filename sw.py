import numpy as np

seq1 = str(input('Choose the first sequence to compare: '))
seq2 = str(input('Choose the second sequence to compare: '))

gap_score = int(input('Set te value for the gap score: '))
match_score = int(input('Set te value for the match score: '))
mismatch_score = int(input('Set te value for the mismatch score: '))


def SW_matrix(seq1, seq2):
    '''returns just an n x m matrix of zeros ''' 
    n = len(seq1)
    m = len(seq2)
    gap = gap_score
    match = match_score
    mismatch = mismatch_score

    matrix = np.zeros((n+1,m+1), dtype=int)

    return matrix

def complete_SWmatrix(): 
    m = SW_matrix(seq1, seq2)
    gap = gap_score
    match = match_score
    mismatch = mismatch_score
    column = 1
    row = 1

    while column < len(seq1) + 1:
        while row < len(m) - 1:
            v1 = m[column][row-1] + gap
            v2 = m[column-1][row] + gap
            if seq1[column-1] == seq2[row-1]:
                v3 = m[column-1][row-1] + match
                values = [v1, v2, v3]
                if max(values) > 0:
                    m[column][row] += max(values)
                else:
                    m[column][row] == 0 #set value to 0 if max(value) < 0
            if seq1[column-1] != seq2[row-1]:
                v3 = m[column-1][row-1] + mismatch
                values = [v1, v2, v3]
                if max(values) > 0:
                    m[column][row] += max(values)
                else:
                    m[column][row] == 0
            row += 1
        row = 1
        column += 1
    return m

def highest():
    '''compute the highest value in the matrix in order to do the local alignment'''
    m = complete_SWmatrix()
    column = 1
    row = 1
    hvalue = [0,[1,1]]

    while column < len(seq1) + 1:
        while row < len(seq2) + 1:
            if m[column][row] > hvalue[0]:
                hvalue[0] = m[column][row]
                hvalue[1] = [column, row]
            row += 1
        row = 1
        column += 1
    return hvalue


def traceback_matrix():
    gap = gap_score
    match = match_score
    mismatch = mismatch_score
    m1 = SW_matrix(seq1, seq2)
    m2 = np.copy(m1)
    column = 1
    row = 1
    
    while column < len(seq1) + 1:
        while row < len(seq2) + 1:
            v1 = m1[column][row-1] + gap
            v2 = m1[column-1][row] + gap
            if seq1[column-1] == seq2[row-1]:
                v3 = m1[column-1][row-1] + match
                values = [v3, v2, v1]
                if max(values) > 0:
                    m1[column][row] += max(values)
                    index = values.index( max(values))
                    if index == 0:
                        m2[column][row] = 1
                    if index == 1:
                        m2[column][row] = 3
                    if index == 2:
                        m2[column][row] = 2
                        
                else:
                    m1[column][row] == 0
                    
            if seq1[column-1] != seq2[row-1]:
                v3 = m1[column-1][row-1] + mismatch
                values = [v3, v2, v1]
                if max(values) > 0:
                    m1[column][row] += max(values)
                    index = values.index( max(values))
                    if index == 0:
                        m2[column][row] = 1
                    if index == 1:
                        m2[column][row] = 3
                    if index == 2:
                        m2[column][row] = 2
                    
                else:
                    m1[column][row] == 0
                    
            row += 1
        row = 1
        column += 1
    return m2

def SW_alignment(seq1, seq2):
    matrix = traceback_matrix()
    hvalue = highest()
    y = hvalue[1][0]
    x = hvalue[1][1]
    aligned1 = ''
    symbols = '' 
    aligned2 = ''
    while y >= 0:
        while x >= 0:
            if matrix[y][x] == 1:
                aligned1 += seq1[y-1]
                symbols += '|'
                aligned2 += seq2[x-1]
                y -= 1
                x -= 1
                
            elif matrix[y][x] == 2:
                aligned1 += '-'
                symbols = ' ' 
                aligned2 += seq2[x-1]
                x -= 1

            elif matrix[y][x] == 3:
                aligned1 += seq1[y-1]
                symbols += ' '
                aligned2 += '-'
                y -= 1
            elif matrix[y][x] <= 0:
                return aligned1[::-1] + '\n'+ symbols[::-1] + '\n' + aligned2[::-1]
        
            


    
