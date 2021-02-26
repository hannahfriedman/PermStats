from permutation import Permutation
import numpy as np
import numpy.linalg as la
import math
from matrixRepCalc import fac

def generate_matrix(n):
    perms = Permutation.group(n)
    permList = []
    i = 0
    nfac = fac(n)
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
    #    print(perm.__str__())
    return mat

def f(n):
    mat = generate_matrix(n)
    print(la.matrix_rank(mat))

