from permutation import Permutation
import matplotlib.pyplot as plt
import math
import seaborn as sns
from typing import Callable
import numpy as np
import perm_stats
import misc
from itertools import combinations
from dft import dft

## This file needs clean up ##

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

def w_ij_kl_unordered_mat(sigma: Permutation, n: int) -> np.array:
    dim = misc.choose(n, 2)
    mat = np.zeros((dim, dim))
    row = 0
    one_n = list(range(1, n+1))
    for i, j in combinations(one_n, 2):
        col = 0
        for k, l in combinations(one_n, 2):
            mat[row, col] = perm_stats.w_ij_kl_unordered(i, j, k, l)(sigma, n)
            col += 1
        row += 1
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

def rho_reg(sigma, n) -> np.array:
    result = np.zeros((math.factorial(n), math.factorial(n)))
    sn = [perm for perm in Permutation.group(n)]
    for row in range(len(sn)):
        for col in range(len(sn)):
            if sn[row] == sigma*sn[col]:
                result[row, col] = 1
                break
    return np.transpose(result)

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
    w_mat = w_ij_kl_mat
    dim = misc.falling_factorial(n, n-2)
    if function == "exced":
        f = perm_stats.excedances
        w_mat = w_ij_mat
        dim = n
        #return representation(perm_stats.excedances, w_ij_mat, n, n, normalize)
    elif function == "fixed":
        f = perm_stats.fixed_points
        w_mat = w_ij_mat
        dim = n
    elif function == "major index":
        f = perm_stats.major_index
        #return representation(perm_stats.major_index, w_ij_kl_mat, misc.falling_factorial(n, n-2), n, normalize)
    elif function == "length":
        f = perm_stats.length
    elif function == "descent":
        f = perm_stats.descent
        #return representation(perm_stats.length, w_ij_kl_mat, misc.falling_factorial(n, n-2), n, normalize)
    else:
        raise("Function not supported.")
    if normalize:
        f = perm_stats.normal(f, n)
    return representation(f, w_mat, dim, n, False)

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
    return representation(perm_stats.w_ij(i,j), w_ij_mat, n, n, False)

def rep_w_ij_kl(i: int, j: int, k: int, l: int, n: int) -> np.array:
    """
    Generates matrix representation of w_ij_kl matrices for convolving
    """
    return representation(perm_stats.w_ij_kl(i,j,k,l), w_ij_kl_mat, misc.falling_factorial(n, n-2), n, False)

# for i in range(4, 7):
#     print(rep_w_ij_kl(2, 3, 4, 1, i))

# for i in range(6, 7):
#     print(np.round(np.real(np.linalg.eig(rep('length', i, False))[0]),3))
    # print(rep('length', i, False))

# print(np.linalg.eig(rep('exced', 5, False))[0])
# print(np.round(np.real(np.linalg.eig(rep('descent', 5, False))[0])))
# print(np.round(np.real(np.linalg.eig(rep('length', 5, False))[0])))
# print(np.round(np.real(np.linalg.eig(rep('major index', 5, False))[0])))

# m = rep('fixed', 5, True)
# factor = rep('fixed', 5, True)
# for _ in range(50):
#     m = m @ factor
# v = np.random.rand(5, 1)
# print(v)
# print(m@v)
# print(m)

# M = np.zeros((5, 5))
# for i in range(1, 6):
#     for j in range(1, 6):
#         M += rep_w_ij(i, j, 5)
# print(M)
    
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

def upper_bound(reps: list) -> float:
    result = 0
    for rep in reps:
        d = rep.shape[0]
        mat = rep @ np.transpose(rep)
        result += d * np.trace(mat)
    return math.sqrt(result/4)

def random_walk(f, k: int, n: int) -> list:
    result = []
    original_reps = dft(perm_stats.normal(f, n),n)[1:]
#    print(original_reps)
    exp_reps = original_reps
    for i in range(1, k+1):
        result.append(upper_bound(exp_reps))
#        print(exp_reps)
        exp_reps = [exp_reps[i] @ original_reps[i] for i in range(len(exp_reps))]
    return result    

def slow_random_walk(f, t: float, k: int, n: int) -> list:
    result = []
    additional = dft(perm_stats.normal(perm_stats.w_ij(1,1), n), n)[1:]
    starting_reps = dft(perm_stats.normal(f, n),n)[1:]
    original_reps = [t*starting_reps[i] + (1-t)*additional[i] for i in range(len(starting_reps))]
    exp_reps = original_reps
    for i in range(1, k+1):
        result.append(upper_bound(exp_reps))
        exp_reps = [exp_reps[i] @ original_reps[i] for i in range(len(exp_reps))]
    return result    


#-----------------------------------------------------------------

def wij_rep(i: int, j: int, n: int) -> np.array:
    i = i-1
    j = j-1
    mat = np.zeros((n,n))
    mat[i, j] = math.factorial(n-1)
    for row in range(mat.shape[0]):
        if row != i:
            for col in range(mat.shape[1]):
                if col != j:
                    mat[row, col] = math.factorial(n-2)
    return mat

print(rep('exced', 4, False))

# n = 3
# mat = sum([wij_rep(i, i, n) for i in range(1, n+1)])/math.factorial(n)

#print(mat)

# n = 5
# k = 20
# y = random_walk(perm_stats.excedances, k, n)
# x = [i for i in range(1, k+1)]
# plt.yscale("linear")
# plt.scatter(x,y)
# plt.title("Exceedences")
# plt.ylabel("Variation Distance")
# plt.xlabel("Number of Steps")
# plt.show()

n = 5
k = 20
y = random_walk(perm_stats.lis, k, n)
x = [i for i in range(1, k+1)]
plt.yscale("log")
plt.scatter(x,y)
plt.title("Longest Increasing Subsequence")
plt.ylabel("Variation Distance")
plt.xlabel("Number of Steps")
plt.savefig('lis_rw.eps')
plt.show()

 
# k = 20
# n = 5
# for t1 in range(1, 20):
#     y = slow_random_walk(perm_stats.length, t1*0.05, k, n)
# #    print(n, np.round_(y, 2))
#     x = [i for i in range(1, k+1)]
#     plt.yscale("linear")
#     plt.scatter(x,y)
#     plt.title("Length")
#     plt.ylabel("Variation Distance")
#     plt.xlabel("Number of Steps")
# plt.show()

# n = 5
# #print(rep("fixed", n, normalize=False))
# mat = rep_w_ij_kl(1, 2, 1, 2, n)
# for i in range(1, n+1):
#     for j in range(i+1, n+1):
#         if not(i == 1 and j == 2):
# #            print(rep_w_ij_kl(i, j, i, j ,n))
#             mat += rep_w_ij_kl(i, j, i, j ,n)
# print(2*mat)
# print(np.round_(np.linalg.eig(2*mat)[0]))
#print(sum([rep_w_ij(i, i, n) for i in range(1, n+1)])
# print(rep_w_ij(1,2,n))
# print(rep("exced", n, normalize=False))
###------------------------------------------------------------------------------------

# Plot excedance random walk with various tau
# for i in range(1,20):
#     pvd(5, "length", 80, 1/i)
# plt.show()

#if __name__ == "__main__":
    # m = rep("exced", 3, False)
    # print(m)
    # ax = sns.heatmap(m, linewidth=0.5)
    # plt.show()
    # count = 0
    # for sigma in Permutation.group(4):
    #     if perm_stats.excedances(sigma, 4) == 2:
    #         count += 1
    # print(count)
