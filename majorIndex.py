from permutation import Permutation
import numpy as np
import numpy.linalg as la
import math

def generate_matrix(n):
    perms = Permutation.group(n)
    permList = []
    i = 0
    nfac = math.factorial(n)
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
    