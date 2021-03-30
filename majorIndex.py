from permutation import Permutation
import numpy as np
import numpy.linalg as la
import math

def fac(n):
    if n == 0:
        return 1
    else:
        return n*fac(n-1)

def generate_matrix(n):
    perms = Permutation.group(n)
    permList = []
    i = 0
    nfac = fac(n)
    mat = np.zeros((n-1, nfac))
    while i < nfac:
        permList.append(next(perms))
        i += 1
    for i in range(1, n):
        for col in range(len(permList)):
            if permList[col](i)>permList[col](i+1):
                    mat[i-1, col] = 1
    print(permList)
    return mat
    
def calc_major_index(perm, n):
    count = 0
    for index in range(1, n):
        if perm(index)>perm(index+1):
            count+= index
    return count

def major_index_tot(n):
    count = 0
    sn = Permutation.group(sn)
    for sigma in sn:
        count += calc_major_index(sigma, n)
    return count



def f(n):
    mat = generate_matrix(n)
    print(la.matrix_rank(mat))