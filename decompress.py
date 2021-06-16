from perm_stats import w_ij
from perm_stats import w_ij_kl
import random_walk
from permutation import Permutation
import numpy as np

from random import *
import math

###------------------------------------------------------------------------------------
def innerProductMatrices(m1: np.matrix, m2: np.matrix) -> float:
    """
    defines the inner product on matrices
    simply multiples the entries of the matrices pointwise, and returns the sum
    matrices act like vectors
    """
    product = 0
    for i in range(m1.shape[0]):
        for j in range(m2.shape[1]):
            product += m1[i, j] * m2[i, j]
    return product

def __main__():
    n = 3
    randMax = 4
    
    snMatrices = {}
    wijMatrices = {}

    for sigma in Permutation.group(n):
        snMatrices[sigma] = random_walk.w_ij_mat(sigma, n)

    linearCombWij = {}
    for i in range(1, n+1):
        for j in range(1, n+1):
            wijMatrices[(i, j)] = random_walk.rep_w_ij(i, j, n) * math.factorial(n-1)
            linearCombWij[(i, j)] = math.floor(random.random()*randMax)
    print("Linear Combo is: \n" + printLinearComboW(linearCombWij, n) + "\n")
    print("Intermediate: \n" + printMatrixSummation(linearCombWij, wijMatrices, n) + "\n")

    resultingMatrix = np.zeros((n, n))
    for pair in linearCombWij.keys():
        resultingMatrix += linearCombWij[pair] * wijMatrices[pair]
    print("Final Result: \n" + str(resultingMatrix) + "\n")

    innerProductDict = {}
    for pair in linearCombWij.keys():
        innerProductDict[pair] = (innerProductMatrices(wijMatrices[pair], resultingMatrix))
    
    modifiedInnerProd = {}
    origMin = min(innerProductDict.values())
    for pair in innerProductDict.keys():
        modifiedInnerProd[pair] = innerProductDict[pair] - origMin

    s = "w_i<-j \n"
    for pair in innerProductDict.keys():
        s += str(pair) + " "
        s += str(linearCombWij[pair]) + " "
        s += str(innerProductDict[pair]) + "\n"
    print(s)

###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
## wij functions

def T_w_ij(n: int) -> np.matrix:
    '''
    Returns an n^2 x n! matrix whose row vectors are wij functions ordered as
    (1,1), (1,2), ..., (1,n), (2, 1), ..., (2, n), ..., (n, 1), ..., (n,n)
    The columns are indexed by elements of Sn
    '''
    result = np.zeros((n**2, math.factorial(n)))

    # we will use these pairs to create our wij functions
    listOfPairs = []
    for i in range(1, n+1):
        for j in range(1, n+1):
            listOfPairs += [(i, j)]
    
    col = 0
    # Loop over columns
    for sigma in Permutation.group(n):
        # Loop over rows
        for row in range(n**2):
            result[row, col] = w_ij(listOfPairs[col][0], listOfPairs[col][1])(sigma, n)
        col += 1

    return result

def TstarT_w_ij(n: int) -> np.matrix:
    '''
    Returns the matrix T*T. Since T is a real matrix, T* is just T transpose.
    '''
    T = Tmatrix_w_ij(n)
    return T.transpose() @ T

def w_ij_vector(i: int, j: int, n: int) -> np.matrix:
    '''
    Returns w_ij as a vector indicator function on Sn
    '''
    M = np.zeros((math.factorial(n), 1))
    row = 0
    w_ij_function = w_ij(i, j)
    # Each entry in the vector represents a permutation
    for sigma in Permutation.group(n):
        M[row, 0] = w_ij_function(sigma, n)
        row += 1
    return M

## Generate T matrices for S3 - S6 ##
# A3 = T_w_ij(3)
# A3t = T_w_ij(3).transpose()
# AtA3 = TstarT_w_ij(3)

# A4 = T_w_ij(4)
# A4t = T_w_ij(4).transpose()
# AtA4 = TstarT_w_ij(4)

# A5 = T_w_ij(5)
# A5t = T_w_ij(5).transpose()
# AtA5 = TstarT_w_ij(5)

# A6 = T_w_ij(6)
# A6t = T_w_ij(6).transpose()
# AtA6 = TstarT_w_ij(6)

def testDecompression(n: int) -> None:
    """
    Test our decompression on NUM_TESTS different random vectors which are linear combinations of wij functions
    """
    Tn = T_w_ij(n)
    TnStar = Tn.transpose()
    a_vector = np.ones((math.factorial(n), 1))
    NUM_TESTS = 30
    
    for runs in range(NUM_TESTS):
        # v is a dictionary with keys being 2-element tuples representing wij functions and keys being the coefficients
        v = {}
        MAX_COEFF = 420
        # Generate a random vector v
        for i in range(1, n+1):
            for j in range(1, n+1):
                v[(i, j)] = randint(0, MAX_COEFF)

        # Convert v into an np array
        v_vector = np.zeros((math.factorial(n), 1))
        for pair in v.keys():
            v_vector += v[pair] * w_ij_vector(pair[0], pair[1], n)

        # Compute the compression u
        u_vector = Tn @ v_vector

        # Decompress u
        sum_u = u_vector.sum()
        prediction = 1/math.factorial(n) * ( (n - 1) * TnStar @ u_vector - (n - 2) * sum_u/n * a_vector )

        # Print statements in case something goes wrong
        if v_vector.all() != prediction.all():
            print("v linear combo is: " + printLinearComboW(v, n))
            print("v vector is: " + np.array2string(v_vector))
            print("prediction is: " + np.array2string(prediction))
            break

    print("All Success!")


###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
## wijkl functions


def T_w_ij_kl(n: int) -> np.matrix:
    '''
    Returns an (n(n-1))^2 x n! matrix whose row vectors are wijkl functions ordered as
    (1,2,1,2), (1,2,1,3), ..., (1,2,1,n), (1,2,2,1), ..., (1,2,n,n-1), (1,3,1,2) ..., (1,n,n,n-1), ..., (2,1,1,2), ..., (n,n-1,n,n-1)
    The columns are indexed by elements of Sn
    '''
    result = np.zeros(((n*(n-1))**2, math.factorial(n)))

     # we will use these 4-element tuples to create our wijkl functions
    listOfQuads = []
    for i in range(1, n+1):
        for j in range(1, n+1):
            if (i != j): # Avoid repetition in first pair
                for k in range(1, n+1):
                    for l in range(1, n+1):
                        if (l != k): # Avoid repetition in second pair
                            listOfQuads += [(i, j, k, l)]
    col = 0
    # Loop over columns
    for sigma in Permutation.group(n):
        # Loop over rows
        for row in range((n*(n-1))**2):
            result[row, col] = w_ij_kl(listOfQuads[row][0], listOfQuads[row][1], listOfQuads[row][2], listOfQuads[row][3])(sigma, n)
        col += 1
    return result

def TstarT_w_ij_kl(n: int) -> np.matrix:
    '''
    Returns the matrix T*T. Since T is a real matrix, T* is just T transpose.
    '''
    T = T_w_ij_kl(n)
    return T.transpose() @ T


def w_ij_kl_vector(i: int, j: int, k:int, l:int, n: int) -> np.matrix:
    '''
    Returns w_ij_kl as a vector indicator function on Sn
    '''
    M = np.zeros((math.factorial(n), 1))
    w_ij_kl_function = w_ij_kl(i, j, k, l)
    # Each entry represents a permutation
    row = 0
    for sigma in Permutation.group(n):
        M[row, 0] = w_ij_kl_function(sigma, n)
        row += 1
    return M

def permutation_representation(sigma: Permutation, n: int) -> np.matrix:
    '''
    Returns the permutation representation (n! x n! matrix) of sigma
    Complexity: O(n!)
    '''
    # Permutations in Sn index the rows and columns of the matrix. We use a dictionary
    # to keep track of which row/column a permutatation is associated with
    sn_dict = {}
    index = 0
    for tau in Permutation.group(n):
        sn_dict[tau] = index
        index += 1
    
    mat = np.zeros((math.factorial(n), math.factorial(n)))
    for tau in Permutation.group(n):
        pi = tau * sigma
        mat[sn_dict[pi], sn_dict[tau]] = 1
    return mat

###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
# Testing our conjecture for wijkl compression for small n

# n = 6
# tstart= TstarT(n)

# Elements that fix n elements (identity)
# I = transposition_matrix(Permutation(), n)

# Elements that fix n-2 elements
# mat2 = np.zeros((math.factorial(n), math.factorial(n)))
# for i in range(1, n):
#     for j in range(i+1,n+1):
#         mat2 += transposition_matrix(Permutation.from_cycles([j,i]), n)

# Elements that fix n-3 elements
# mat3 = np.zeros((math.factorial(n), math.factorial(n)))
# for i in range(1, n-1):
#     for j in range(i+1,n):
#         for k in range(j+1, n+1):
#             mat3 += transposition_matrix(Permutation.from_cycles([i,j,k]), n)
#             mat3 += transposition_matrix(Permutation.from_cycles([i,k,j]), n)

# Elements that fix n-4 elements
# mat4 = np.zeros((math.factorial(n), math.factorial(n)))
# for i in range(1, n-2):
#     for j in range(i+1,n-1):
#         for k in range(j+1, n):
#             for l in range(k+1, n+1):
#                 mat4 += transposition_matrix(Permutation.from_cycles([i,j,k,l]), n)
#                 mat4 += transposition_matrix(Permutation.from_cycles([i,j,l,k]), n)
#                 mat4 += transposition_matrix(Permutation.from_cycles([i,k,j,l]), n)
#                 mat4 += transposition_matrix(Permutation.from_cycles([i,k,l,j]), n)
#                 mat4 += transposition_matrix(Permutation.from_cycles([i,l,j,k]), n)
#                 mat4 += transposition_matrix(Permutation.from_cycles([i,l,k,j]), n)

#                 mat4 += transposition_matrix(Permutation.from_cycles([i,j], [k,l]), n)
#                 mat4 += transposition_matrix(Permutation.from_cycles([i,k], [j,l]), n)
#                 mat4 += transposition_matrix(Permutation.from_cycles([i,l], [j,k]), n)

# print(np.array_equal(tstart - 2*mat2 - 6*I, np.zeros((math.factorial(n), math.factorial(n)))))
# print(np.array_equal(tstart - 6*mat2 - 2*mat3 - 12*I, np.zeros((math.factorial(n), math.factorial(n)))))
# print(np.array_equal(tstart - 12*mat2 - 6*mat3 - 2*mat4 - 30*I, np.zeros((math.factorial(n), math.factorial(n)))))

###------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------
## Print Functions

def printOneLineNotation(sigma: Permutation, n: int) -> str:
    """
    prints the permutation sigma using one line notation
    ex: [3 4 5 1 2]
    """
    s = "["
    for i in range(1, n+1):
        s += str(sigma(i))
        s += " "
    s = s[:-1] # remove last space
    s += "]"
    return s


def printLinearComboW(d: dict, n: int) -> str:
    """
    takes a of w_ij's
    the dictionary stores keys that are w_ij's, and the values are the coefficient
    prints it nicely
    """
    s = ""
    for pair in d.keys():
        s += str(d[pair])
        s += "(w_" + str(pair[0]) + "<-" + str(pair[1]) + ")"
        s += " + "
    s = s[:-3] # remove last plus
    return s

def printMatrixSummation(linearCombo: dict, matrixDict: dict, n: int) -> str:
    """
    takes a linear combination of matrices, and their w_ij's
    both linearCombo and matrixDict use w_ij's as keys
         the linearCombo associates each key with its coefficient
         the matrixDict associates each key with its representation matrix
    prints the linear combination nicely
    """
    s = ""
    for pair in linearCombo.keys():
        s += str(linearCombo[pair])
        s += str(matrixDict[pair])
        s += " + \n"
    s = s[:-4] # removes last plus
    return s
###------------------------------------------------------------------------------------
