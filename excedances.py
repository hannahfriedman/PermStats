from permutation import Permutation
import numpy as np
import numpy.linalg as la
import math
from wij import w_ij

def generate_matrix(n):
    perms = Permutation.group(n)
    permList = []
    i = 0
    nfac = math.factorial(n)
    mat = np.zeros((n-1, nfac))
    while i < nfac:
        permList.append(next(perms))
        i += 1
    for row in range(n-1):
        for col in range(len(permList)):
            for i in range(1, n+1):
                if permList[col].__call__(i)==row + 2  and row + 2>i:
                    mat[row, col] = 1
    #for perm in permList:
    #   print(perm.__str__())
    return mat

def sn_act_on_matrix(n, i):
    """
    returns the matrix of excedance indicator functions for Sn after being acted on by (i,i+1)
    """
    mat = generate_matrix(n)
    nfac = math.factorial(n)
    perms = Permutation.group(n)
    switched = []
    #These dictionaries allow us to access permutations by the index they represent and vice versa
    perm_index_dict = {}
    index_perm_dict = {}
    perm_list = []
    #Populate the dictionaries
    for j in range(nfac):
        #The following line is printing and i'm not sure how to make it stop
        p = next(perms)
        perm_index_dict[p] = j
        index_perm_dict[j] = p
        perm_list.append(p)
    for col in range(len(mat[0])):
        perm = index_perm_dict[col]
        if perm not in switched:
            perm_acted_on = perm*Permutation.cycle(i, i+1)
            col2 = perm_index_dict[perm_acted_on]      
            perm_list[col], perm_list[col2] = perm_list[col2], perm_list[col]
            for row in range(n-1):
                mat[(row, col)], mat[(row, col2)] = mat[(row, col2)], mat[(row, col)]
            switched.append(perm_acted_on)
    #for perm in perm_list:
    #   print(perm.__str__())
    return mat
    
def mega_matrix(n):
    """
    what does this do?
    """
    nfac = math.factorial(n)
    mat = np.zeros((n**2 - n, nfac))
    for matrix_index in range(n):
        if matrix_index == 0:
            matrix = generate_matrix(n)
        else:
            matrix = sn_act_on_matrix(n, matrix_index)
        for row in range(n-1):
                mat[row + (n-1)*matrix_index] = matrix[row]
    return mat

