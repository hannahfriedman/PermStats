from excedances import count_excedances
from excedances import total_exced
from permutation import Permutation
import numpy as np
from wij import w_ij
from misc import matlab_syntax

def w_ij_mat(n, sigma):
    mat = np.zeros((n,n))
    for i in range(1, n+1):
        for j in range(1, n+1):
            mat[i-1, j-1] = w_ij(sigma, j, i)
    return mat

def w_ij_mat_sn(n, sn):
    mats = {}
    for sigma in sn:
        mats[sigma] = w_ij_mat(n, sigma)
    return mats

def representation(n, function, function_sum):
    sn = Permutation.group(n)
    mats = w_ij_mat_sn(n, sn)
    tot = function_sum(n)
    mat = np.zeros((n,n))
    for sigma in sn:
        factor = function(sigma)/tot
        print(factor * mats[sigma])
        mat = mat + factor * mats[sigma]
    return mat

print(representation(4, count_excedances, total_exced))


