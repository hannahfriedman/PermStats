from permutation import Permutation
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Callable

import numpy as np
import perm_stats
import misc

###------------------------------------------------------------------------------------
#Generate w_ij, w_ij_kl matrices for Sn
def w_ij_mat(sigma: Permutation, n: int) -> np.array:
    '''
    Create an n x n representation of sigma
    This representation contains the trivial and (n-1, 1) representations
    '''
    mat = np.zeros((n,n))
    for i in range(1, n+1):
        for j in range(1, n+1):
            # i represents the rows, j represents the columns
            mat[i-1, j-1] = perm_stats.w_ij(i,j)(sigma,n)
    return mat

def w_ij_kl_mat(sigma: Permutation, n: int) -> np.array:   
    '''
    Create a representation of sigma
    This representation contains the (n), (n-1, 1), (n-2, 2), (n-2, 1, 1) representations
    '''
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

def w_mat_sn(w: Callable[[Permutation, int], np.array], n: int) -> dict:
    """
    Returns a dictionary of w-representations of permutations in Sn 
    w should be either w_ij_mat or w_ij_kl_mat
    """
    sn = Permutation.group(n)
    mats = {}
    for sigma in sn:
        mats[sigma] = w(sigma, n)
    return mats
###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
# Generate matrix representations of functions on Sn
def representation(function: Callable[[Permutation, int], int], 
                   w: Callable[[Permutation, int], np.array], 
                   dim: int, 
                   n: int, 
                   normalize: bool) -> np.array:
    """
    Returns a linear combination of w_ij or w_ij_kl matrices 
    with coefficients being function of the permutation represented
    by the matrix
    function is a permutation statistic
    w indicates whether to use the w_ij or w_ij_kl representation
    dim is the dimension of the matrices, should be n for w_ij and n falling factorial 2 for w_ij_kl
    n is the n in Sn 
    """
    sn = Permutation.group(n)
    # Create a dictionary mapping permutations to their w-matrices
    mats = w_mat_sn(w, n)
    # Sum function evaluated for every permutation--used to scale so all rows
    # and columns of matrices sum to 1
    if normalize:
        tot = perm_stats.total(function, n)
    else:
        tot = 1
    mat = np.zeros((dim,dim))
    # Add matrices to mat after scaling them appropriately 
    for sigma in sn:
        factor = function(sigma, n)/tot
        mat = mat + factor * mats[sigma]
    return mat

def rep(function: str, n: int, normalize=True) -> np.array:
    '''
    User friendly version of representation function
    '''
    if function == "exced":
        return representation(perm_stats.excedances, w_ij_mat, n, n, normalize)
    if function == "major index":
        return representation(perm_stats.major_index, w_ij_kl_mat, misc.falling_factorial(n, n-2), n, normalize)
    if function == "length":
        return representation(perm_stats.length, w_ij_kl_mat, misc.falling_factorial(n, n-2), n, normalize)

def rep_tau(tau: float, function: Callable[[Permutation, int], int], n: int) -> np.array:
    """
    Returns a "watered-down" representation of function, where tau times the representation of the 
    function is added to (1-tau) times the representation of w_11 or w_12_12
    tau       float between 0 and 1, inclusive
    function  permutation statistic
    n         n in Sn
    """
    # Create representation of function scaled by tau
    mat1 = tau*rep(function, n)
    # Create representation of w-matrix
    if function == "exced":  # Use w_ij if function is excedances...
        w = rep_w_ij(1,1, n)
    else:                    # ...use w_ij_kl otherwise
        w = rep_w_ij_kl(1,2,1,2,n)
    # Scale representation of w-matrix by 1-tau
    mat2 = (1-tau)*w
    return mat1+mat2

def rep_w_ij(i: int, j: int, n: int) -> np.array:
    """
    Generates matrix representation of w_ij matrices for convolving
    """
    return representation(perm_stats.w_ij(i,j), w_ij_mat, n, n)

def rep_w_ij_kl(i: int, j: int, k: int, l: int, n: int) -> np.array:
    """
    Generates matrix representation of w_ij_kl matrices for convolving
    """
    return representation(perm_stats.w_ij_kl(i,j,k,l), w_ij_kl_mat, misc.falling_factorial(n, n-2), n)
###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
# Variation distance and Exponentiating for Random Walks
def var_dist(mat: np.array) -> float:
    """
    computes variation distance of a matrix 
    the variation distance is the maximum difference between the value in the ith row and jth column
    and the norm_dist_val, which is 1 over the dimension of the matrix
    We use this value to compare because a normal distribution means every w_ij matrix gets chosen with 
    equal probability, so we have the same value in each entry. However, we still need each row and column to
    sum to 0, so the norm_dist_val is 1 over the dimension of the matrix.
    """
    norm_dist_val = 1/mat.shape[0] 
    max_dif = 0
    for row in range(0, mat.shape[0]):
        for col in range(0, mat.shape[1]):
            # If the difference between an entry and norm_dist_val is greater than the maximum difference...
            if abs(mat[row, col] - norm_dist_val) > max_dif:
                # Set maximum difference to be that difference
                max_dif = abs(mat[row, col] - norm_dist_val)
    return max_dif

def plot_var_dist(mat: np.array, max_pow: int) -> list:
    """
    Plots variation distance of an array when raised to powers from 1 to max_pow
    """
    powers = []
    # Fill powers with variation distances of different powers
    for power in range(1, max_pow+1):
        powers.append(var_dist(np.linalg.matrix_power(mat, power)))
    # Plot powers using pyplot
    plt.scatter(range(1, max_pow+1), powers)
    plt.yscale("linear")
    return powers

def pvd(function: Callable[[Permutation, int], int], 
        max_pow: int, 
        tau: float, 
        n: int) -> None:
    """
    User friendly version of plot_var_dist
    function - permutation statistic
    max_pow  - greatest power to which to raise matrix
    tau      - parameter to control how quickly variation distance goes to 0
               use tau = 1 for the unfiltered permutation statistic
    n        - n in Sn
    """
    mat = rep_tau(tau, function, n)
    plot_var_dist(mat, max_pow)

def pvd_w_ij(i: int, j: int, max_pow: int, n: int) -> None:
    """
    User friendly version of plot_var_dist for w_ij functions
    """
    mat1 = rep_w_ij(i, j, n)
    plot_var_dist(mat1, max_pow)

def pvd_w_ij_kl(i: int, j: int, k: int, l: int, max_pow: int, n: int) -> None:
    """
    User friendly version of plot_var_dist for w_ij functions
    """
    mat1 = rep_w_ij(i, j, n)/2 + rep_w_ij(k,l, n)/2
    plot_var_dist(mat1, max_pow)
    
###------------------------------------------------------------------------------------

# Plot excedance random walk with various tau
# for i in range(1,20):
#     pvd(5, "length", 80, 1/i)
# plt.show()

if __name__ == "__main__":
    m = rep("length", 4, False)
    print(m)
    # ax = sns.heatmap(m, linewidth=0.5)
    # plt.show()