from matrixRepCalc import matrix_rep
from matrixRepCalc import fac
from matrixRepCalc import adjust_zeros
from scipy.linalg import block_diag
from excedances import count_excedances
from excedances import total_exced
from misc import matlab_syntax
import numpy as np
import math
from permutation import Permutation


#To do: we care about the matrix phi* m^n phi





def normalized_DFT_excedances(n):
    """
    n--int n in S_n
    """
    nfac = fac(n)
    sn = Permutation.group(n)
    dft = []
    rho = matrix_rep(n)
    for i in range(nfac):
        perm = next(sn)
        scaled_excedances = count_excedances(perm, n)/total_exced(n)
        if i == 0:
            for mat in rho[perm]:
                dft.append(scaled_excedances*mat)
        else:
            for i in range(len(rho[perm])):
                dft[i] = np.add(dft[i], scaled_excedances*rho[perm][i])
    return adjust_zeros(dft)


def generate_Mk_exced(matk):
    return block_diag(*([matk] * len(matk)))

def generate_M_exced(DFT):
    m_k = [generate_Mk_exced(matk) for matk in DFT]
    return block_diag(*m_k)



def psi_k_s(n, matk):
    """
    Takes in a matrix in the matrix representation of s
    """
    d_k = len(matk)
    cardinality = fac(n)
    col_vec = np.zeros((d_k**2))
    scale_fac = math.sqrt(d_k/cardinality)
    for col in range(d_k):
        for row in range(d_k):
            col_vec[(row-1)*d_k + col] = scale_fac*matk[row, col]
    return col_vec

def phi_s(n, mat_rep):
    psi_s = [psi_k_s(n, matk) for matk in mat_rep]
    return np.concatenate((psi_s))

def phi(n):
    phis = [phi_s(n, mat_rep) for mat_rep in matrix_rep(n).values()]
    return np.array(phis)

def transition_matrix_exced(n):
    dft = normalized_DFT_excedances(n)
    m = generate_M_exced(dft)
    phi_ = phi(n)
    phi_star = np.linalg.inv(phi_)
    print(m)
    return adjust_zeros([np.matmul(phi_star, np.matmul(m, phi_))])[0]


def __main__():
    trans_mat = transition_matrix_exced(3)
    trans_mat1 = transition_matrix_exced(4)
    #trans_mat2 = transition_matrix_exced(5)
    #print(trans_mat)
    #print(trans_mat1)
    #print(trans_mat2)
    #print(matlab_syntax(trans_mat2))

__main__()