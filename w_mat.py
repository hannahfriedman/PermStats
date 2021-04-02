from excedances import count_excedances
from excedances import total_exced
from majorIndex import calc_major_index
from permutation import Permutation
from wij import w_ij
from wij import w_ij_kl
from misc import matlab_syntax
import matplotlib.pyplot as plt
import numpy as np

###------------------------------------------------------------------------------------
#Permutation statistic function I can't impor for some reason

def total_length(n):
    """
    returns the total length (number of inversion)
    over all permutations in Sn
    """
    sn = Permutation.group(n)
    length = 0
    for sigma in sn:
        length += sigma.inversions()
    return length

def length(sigma, n):
    return sigma.inversions()

def major_index_tot(n):
    """
    returns the total major index
    over all permutations in Sn
    """
    count = 0
    sn = Permutation.group(n)
    for sigma in sn:
        count += calc_major_index(sigma, n)
    return count
###------------------------------------------------------------------------------------

###------------------------------------------------------------------------------------
#Generate w_ij, w_ij_kl matrices for Sn
def w_ij_mat(n, sigma):
    mat = np.zeros((n,n))
    for i in range(1, n+1):
        for j in range(1, n+1):
            mat[i-1, j-1] = w_ij(sigma, j, i)
    return mat

def w_ij_kl_mat(n, sigma):
    mat = np.zeros((n**2,n**2))
    for i in range(1, n+1):
        for j in range(1, n+1):
            for k in range(1, n+1):
                for l in range(1, n+1):
                    row = (i-1)*n + (j-1)
                    col = (k-1)*n + (l-1)
                    mat[row, col] = w_ij_kl(sigma, i, j, k, l)
    return mat

def w_mat_sn(n, w):
    sn = Permutation.group(n)
    mats = {}
    for sigma in sn:
        mats[sigma] = w(n, sigma)
    return mats
###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
# Linear combinations of w matrices using coefficients
def representation(n, function, function_sum, w, dim):
    sn = Permutation.group(n)
    mats = w_mat_sn(n, w)
    tot = function_sum(n)
    mat = np.zeros((dim,dim))
    for sigma in sn:
        factor = function(sigma, n)/tot
        mat = mat + factor * mats[sigma]
    return mat

def rep(n, function):
    '''
    User friendly version of representation function
    '''
    if function == "exced":
        return representation(n, count_excedances, total_exced, w_ij_mat, n)
    if function == "major index":
        return representation(n, calc_major_index, major_index_tot, w_ij_kl_mat, n**2)
    if function == "length":
        return representation(n, length, total_length, w_ij_kl_mat, n**2)
###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
# Convert the matrices into complex, diagonal matrices and analyze
def diagonalize(n, rep):
    mat = np.zeros((n,n), dtype = complex)
    evals = np.linalg.eig(rep)[0]
    for i in range(n):
        mat[i,i] = evals[i]
    return mat

def var_dist(mat):
    max_val = -1
    for row in range(1, mat.shape[0]):
        for col in range(1, mat.shape[1]):
            if norm(mat[row,col]) > max_val:
                max_val = norm(mat[row,col])
    return max_val

def norm(c):
    return c.real**2 + c.imag**2

def plot_var_dist(mat, max_pow):
    powers = []
    for power in range(1, max_pow+1):
        powers.append(var_dist(np.linalg.matrix_power(mat, power)))
    plt.yscale('log')
    plt.scatter(range(1, max_pow+1), powers)
    plt.show()
    return powers

def pvd(n, function, max_pow):
    """
    User friendly version of plot_var_dist
    """
    mat = rep(n, function)
    if function == "exced":
        dim = n
    else:
        dim = n**2
    print(matlab_syntax(mat))
    diag = diagonalize(dim, mat)
    #print(matlab_syntax(diag))
    plot_var_dist(diag, max_pow)

M = rep(4, "exced")
N = rep(4, "major index")
O = rep(4, "length")
#print(M)
pvd(4, "length", 20)
# print(O)



# plot_var_dist(diag, 20)

# evals, evecs = np.linalg.eig(rep)
# print(evals)
# print(evecs)


#print(matlab_syntax(rep))









