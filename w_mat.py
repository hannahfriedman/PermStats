from permutation import Permutation
import matplotlib.pyplot as plt
import numpy as np
import perm_stats
import misc

###------------------------------------------------------------------------------------
#Generate w_ij, w_ij_kl matrices for Sn
def w_ij_mat(sigma: Permutation, n: int) -> np.array:
    mat = np.zeros((n,n))
    for i in range(1, n+1):
        for j in range(1, n+1):
            mat[i-1, j-1] = perm_stats.w_ij(i,j)(sigma,n)
    return mat

def w_ij_kl_mat(sigma: Permutation, n: int) -> np.array:   
    dim = misc.falling_factorial(n, n-2) 
    mat = np.zeros((dim, dim))
    shift_ij = 1  # keep track of how much our index differs if we used (i,i) pairs
    for i in range(1, n+1):
        for j in range(1, n+1):
            row = (i-1)*n + j -  shift_ij # compute the correct row as we would if we included (i,i) pairs and subtract the shift
            shift_kl = 1
            if i != j:
                for k in range(1, n+1):
                    for l in range(1, n+1):
                        if k != l:
                            col = (k-1)*n + l -  shift_kl
                            mat[row, col] = perm_stats.w_ij_kl(i, j, k, l)(sigma, n)
                        else:
                            shift_kl += 1
            else:
                shift_ij += 1
    return mat

def w_mat_sn(n, w):
    sn = Permutation.group(n)
    mats = {}
    for sigma in sn:
        mats[sigma] = w(sigma, n)
    return mats
###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
# Linear combinations of w matrices using coefficients
def representation(n, function, w, dim):
    sn = Permutation.group(n)
    mats = w_mat_sn(n, w)
    tot = perm_stats.total(function, n)
    mat = np.zeros((dim,dim))
    for sigma in sn:
        factor = function(sigma, n)/tot
        mat = mat + factor * mats[sigma]
    return mat

def rep(n, function):
    '''
    #User friendly version of representation function
    '''
    if function == "exced":
        return representation(n, perm_stats.excedances, w_ij_mat, n)
    if function == "major index":
        return representation(n, perm_stats.major_index, w_ij_kl_mat, misc.falling_factorial(n, n-2))
    if function == "length":
        return representation(n, perm_stats.length, w_ij_kl_mat, misc.falling_factorial(n, n-2))

def rep_w_ij(i: int, j: int, n: int) -> np.array:
    """
    Generates linear combinations of w_ij matrices for convolving
    """
    return representation(n, perm_stats.w_ij(i,j), w_ij_mat, n)

###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------

def convolve(n: int) -> list:
    mats = []
    for i in range(1, n+1):
        for j in range(1, n+1):
            ij_mat = rep_w_ij(i, j, n)
            for k in range(1, n+1):
                for l in range(1, n+1):
                    kl_mat = rep_w_ij(k, l, n)
                    print("w" + str(i) + str(j) + " * w" + str(k) + str(l) + ":")
                    print(np.matmul(ij_mat, kl_mat))
                    mats.append(np.matmul(ij_mat, kl_mat))
    return mats

def count_unique(mats: list) -> list:
    used = []
    for mat in mats:
        in_list = False
        for x in used:
            if (mat == x).all():
                in_list = True
                break
        if not in_list:
            used.append(mat)

    return used


def print_mult_table(n: int) -> None:
    flat = convolve(n)
    arr = []
    for i in range(n**2):
        arr.append(flat[i*(n**2):(i+1)*(n**2)])
    for x in arr:
        for y in x:
            print(y)
        print('\n')



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
        dim = misc.falling_factorial(n, n-2)
    print(misc.matlab_syntax(mat))
    diag = diagonalize(dim, mat)
    #print(misc.matlab_syntax(diag))
    plot_var_dist(diag, max_pow)

M = rep(4, "exced")
N = rep(3, "major index")
O = rep(4, "length")













