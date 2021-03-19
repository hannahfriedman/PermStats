from matrixRepCalc import DFT_excedances
from matrixRepCalc import matrix_rep
from matrixRepCalc import fac
from scipy.linalg import block_diag
import numpy as np
import math
from permutation import Permutation


#To do: Replace DFT is normalized function and implement transition matrix function


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




