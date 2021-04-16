from perm_stats import w_ij
import w_mat
from permutation import Permutation
import numpy as np

import random
import math

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
        snMatrices[sigma] = w_mat.w_ij_mat(sigma, n)

    linearCombWij = {}
    for i in range(1, n+1):
        for j in range(1, n+1):
            wijMatrices[(i, j)] = w_mat.rep_w_ij(i, j, n) * math.factorial(n-1)
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

def Tmatrix(n: int) -> np.matrix:
    result = np.zeros((math.factorial(n), n**2))
    listOfPairs = []
    for i in range(1, n+1):
        for j in range(1, n+1):
            listOfPairs += [(i, j)]
    
    row = 0
    for sigma in Permutation.group(n):
        for col in range(n**2):
            result[row, col] = w_ij(listOfPairs[col][0], listOfPairs[col][1])(sigma, n)
        row += 1

    return result

def TadjointTmat(n: int) -> np.matrix:
    T = Tmatrix(n)
    return np.matmul(T.transpose(), T)
