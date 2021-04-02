from permutation import Permutation
import numpy as np
from excedances import count_excedances
from majorIndex import calc_major_index
from misc import matlab_syntax

# Excedance Transition Matrix S3 has eigenvaluse 6 with eigenvector all 1s, rank = 5


def s_n_list(n):
    return [*Permutation.group(n)]

def rep_k(n, sigma, s_n):
    nfac = len(s_n)
    mat = np.zeros((nfac, nfac))
    for tau_index in range(nfac):
        sigma_tau = sigma * s_n[tau_index]
        sigma_tau_index = s_n.index(sigma_tau)
        mat[tau_index, sigma_tau_index] = 1
    return mat

def regular_representation(n):
    s_n = s_n_list(n)
    rho = {}
    for sigma in s_n:
        rho[sigma] = rep_k(n, sigma, s_n)
    return rho

def transition_excedances(n):
    s_n = s_n_list(n)
    reg_rep = regular_representation(n)
    trans_mat = np.zeros((len(s_n), len(s_n)))
    excedances_total = 0
    for sigma in s_n:
        excedances = count_excedances(sigma, n)
        trans_mat = trans_mat + excedances*reg_rep[sigma]
        excedances_total += excedances
    return [trans_mat, excedances_total]

def transition_major_index(n):
    s_n = s_n_list(n)
    reg_rep = regular_representation(n)
    trans_mat = np.zeros((len(s_n), len(s_n)))
    major_index_total = 0
    for sigma in s_n:
        major_index = calc_major_index(sigma, n)
        trans_mat = trans_mat + major_index*reg_rep[sigma]
        major_index_total += major_index
    return [trans_mat, major_index_total]





    
    
def __main__():
    trans_mat = transition_major_index(3)[0]
    #print(transition_major_index(5))
    #print(wolfram_syntax(trans_mat))
    print(matlab_syntax(trans_mat))
    #print(transition_excedances(4))


    

__main__()