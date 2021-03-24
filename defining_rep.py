from permutation import Permutation
import numpy as np
import math
from matrixRepCalc import adjust_zeros


def rho(n):
    sn = Permutation.group(n)
    d = {}
    for sigma in sn:
        d[sigma] = gen_rep(sigma, n)
    return d

def gen_rep(sigma, n):
    mat = np.zeros((n,n))
    for row in range(0, n):
        mat[row, sigma(row+1)-1] = 1
    return mat

def cob_s_three():
    v1 = [1/math.sqrt(3),1/math.sqrt(3),1/math.sqrt(3)]
    v2 = [-1/math.sqrt(2),0,1/math.sqrt(2)]
    v3 = [1/math.sqrt(6),-2/math.sqrt(6),1/math.sqrt(6)]
    return np.array([v1,v2,v3])

def print_images():
    sn = Permutation.group(3)
    p = rho(3)
    cob = cob_s_three()
    for sigma in sn:
        print(adjust_zeros([np.matmul(cob, np.matmul(p[sigma], np.transpose(cob)))])[0])
