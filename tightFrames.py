from perm_stats import w_ij
from perm_stats import w_ij_kl
import w_mat
from permutation import Permutation
import numpy as np

from random import *
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

def TmatTstar(n: int) -> np.matrix:
    T = Tmatrix(n)
    return np.matmul(T, T.transpose())

def w_ij_vector(i: int, j: int, n: int) -> np.matrix:
    M = np.zeros((math.factorial(n), 1))
    row = 0
    w_ij_function = w_ij(i, j)
    for sigma in Permutation.group(n):
        M[row, 0] = w_ij_function(sigma, n)
        row += 1
    return M

A3 = Tmatrix(3).transpose()
A3t = Tmatrix(3)
AtA3 = np.matmul(A3t, A3)

A4 = Tmatrix(4).transpose()
A4t = Tmatrix(4)
AtA4 = np.matmul(A4t, A4)

A5 = Tmatrix(5).transpose()
A5t = Tmatrix(5)
AtA5 = np.matmul(A5t, A5)

A6 = Tmatrix(6).transpose()
A6t = Tmatrix(6)
AtA6 = np.matmul(A6t, A6)

def testDecompression(n: int) -> None:
    """
    testing our decompression formula to make sure it works
    """
    
    TnStar = Tmatrix(n)
    Tn = TnStar.transpose()
    a_vector = np.ones((math.factorial(n), 1))
    
    for runs in range(30):
        v = {}
        MAX_COEFF = 420
        for i in range(1, n+1):
            for j in range(1, n+1):
                v[(i, j)] = randint(0, MAX_COEFF)

        v_vector = np.zeros((math.factorial(n), 1))
        for pair in v.keys():
            v_vector += v[pair] * w_ij_vector(pair[0], pair[1], n)

        u_vector = Tn @ v_vector

        sum_u = u_vector.sum()

        prediction = 1/math.factorial(n) * ( (n - 1) * TnStar @ u_vector - (n - 2) * sum_u/n * a_vector )

        # print("v linear combo is: " + printLinearComboW(v, n))

        if v_vector.all() != prediction.all():
            print("v linear combo is: " + printLinearComboW(v, n))
            print("v vector is: " + np.array2string(v_vector))
            print("prediction is: " + np.array2string(prediction))
            break

    print("All Success!")



def TmatrixW_ij_kl(n: int) -> np.matrix:
    result = np.zeros(((n*(n-1))**2, math.factorial(n)))
    listOfQuads = []
    for i in range(1, n+1):
        for j in range(1, n+1):
            if (i != j):
                for k in range(1, n+1):
                    for l in range(1, n+1):
                        if (l != k):
                            listOfQuads += [(i, j, k, l)]
    col = 0
    for sigma in Permutation.group(n):
        for row in range((n*(n-1))**2):
            result[row, col] = w_ij_kl(listOfQuads[row][0], listOfQuads[row][1], listOfQuads[row][2], listOfQuads[row][3])(sigma, n)
        col += 1
    return result

def TstarT(n: int) -> np.matrix:
    T = TmatrixW_ij_kl(n)
    return T.transpose() @ T

# for eig in np.linalg.eigh(TstarT(6))[0]:
    # print(round(eig))

# for eig in np.linalg.eigh(TmatTstar(6))[0]:
    # print(round(eig))

def uniqueEigs(m: np.matrix) -> list:
    dictOfEigs = {}
    for eig in np.linalg.eigh(m)[0]:
        roundedEig = round(eig)
        if roundedEig not in dictOfEigs.keys():
            dictOfEigs[roundedEig] = 1
        else:
            dictOfEigs[roundedEig] += 1
    return dictOfEigs



def w_ij_kl_vector(i: int, j: int, k:int, l:int, n: int) -> np.matrix:
    M = np.zeros((math.factorial(n), 1))
    row = 0
    w_ij_kl_function = w_ij_kl(i, j, k, l)
    for sigma in Permutation.group(n):
        M[row, 0] = w_ij_kl_function(sigma, n)
        row += 1
    return M

vec1 = w_ij_kl_vector(1,2,3,4,5)
vec2 = w_ij_kl_vector(3,4,1,2,4)
vec3 = w_ij_kl_vector(3,2,3,1,4)

tstart_four = TstarT(5)


print(tstart_four @ vec1)
# print(tstart_four @ vec2)
# print(tstart_four @ vec3)