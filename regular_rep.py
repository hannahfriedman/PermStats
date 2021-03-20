from permutation import Permutation
import numpy as np
from excedances import count_excedances

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

def wolfram_syntax(mat):
    string = "{"
    for row in range(len(mat)):
        string+="{"
        for col in range(len(mat[0])-1):
            string += str(mat[row, col])
            string += ","
        string += str(mat[row, len(mat[0])-1])
        if row == len(mat) - 1:
            string += "}"
        else:
            string += "},"
    string+= "}"
    print(string)

def matlab_syntax(mat):
    string = "["
    for row in range(len(mat)):
        for col in range(len(mat)):
            string += " "
            string += str(mat[row, col])
        if row == len(mat) - 1:
            string += "]"
        else:
            string += ";"
    print(string)



    

    
def __main__():
    trans_mat = transition_excedances(6)[0]
    #print(wolfram_syntax(trans_mat))
    print(matlab_syntax(trans_mat))
    #print(transition_excedances(4))


    

__main__()