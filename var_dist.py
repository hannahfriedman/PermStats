from matrixRepCalc import DFT_excedances
from matrixRepCalc import DFT_major_index
from permutation import Permutation
from excedances import count_excedances
import numpy as np


def total_exced(n):
    sn = Permutation.group(n)
    count = 0
    for sigma in sn:
        count += count_excedances(sigma, n)
    return count

def normalize(mat, stat):
    for row in range(mat.shape[0]):
        for col in range(mat.shape[1]):
            mat[row, col] =  mat[row, col]/stat
    return mat

def var_dist(mat):
    max_val = -1
    for row in range(mat.shape[0]):
        for col in range(mat.shape[1]):
            if abs(mat[row,col]) > max_val:
                max_val = abs(mat[row,col])
    return max_val

def __main__():
    mat = DFT_excedances(5)[2]
    mat = normalize(mat, total_exced(5))
    print(var_dist(mat))
    mat2 = np.linalg.matrix_power(mat, 2)
    print(var_dist(mat2))
    mat3 = np.linalg.matrix_power(mat, 3)
    print(var_dist(mat3))
    mat4 = np.linalg.matrix_power(mat, 4)
    print(var_dist(mat4))
    mat5 = np.linalg.matrix_power(mat, 5)
    print(var_dist(mat5))
    mat6 = np.linalg.matrix_power(mat, 6)
    print(var_dist(mat6))
    mat7 = np.linalg.matrix_power(mat, 7)
    print(var_dist(mat7))
    mat8 = np.linalg.matrix_power(mat, 8)
    print(var_dist(mat8))
    

    '''
    mat = DFT_major_index(5)[2]
    mat = normalize(mat, major_index(5))
    print(mat)
    print(var_dist(mat))
    mat2 = np.linalg.matrix_power(mat, 2)
    print(var_dist(mat2))
    mat3 = np.linalg.matrix_power(mat, 3)
    print(var_dist(mat3))
    mat4 = np.linalg.matrix_power(mat, 4)
    print(var_dist(mat4))
    mat5 = np.linalg.matrix_power(mat, 5)
    print(var_dist(mat5))
    mat6 = np.linalg.matrix_power(mat, 6)
    print(var_dist(mat6))
    mat7 = np.linalg.matrix_power(mat, 7)
    print(var_dist(mat7))
    mat8 = np.linalg.matrix_power(mat, 8)
    print(var_dist(mat8))
    '''


__main__()

