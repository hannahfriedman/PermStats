from matrixRepCalc import DFT_excedances
from matrixRepCalc import DFT_major_index
from matrixRepCalc import DFT_length
from matrixRepCalc import DFT_w_ij
from permutation import Permutation
from excedances import count_excedances
from misc import total_length
import matplotlib.pyplot as plt
import numpy as np


def var_dist(mat):
    max_val = -1
    for row in range(mat.shape[0]):
        for col in range(mat.shape[1]):
            if abs(mat[row,col]) > max_val:
                max_val = abs(mat[row,col])
    return max_val

def plot_var_dist(mat, max_pow):
    powers = []
    for power in range(1, max_pow+1):
        powers.append(var_dist(np.linalg.matrix_power(mat, power)))
    plt.yscale('log')
    plt.scatter(range(1, max_pow+1), powers)
    plt.show()
    return powers

def plot_var_dist2(mat, mat1, max_pow):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.yscale('log')
    powers = []
    powers1 = []
    for power in range(1, max_pow+1):
        powers.append(var_dist(np.linalg.matrix_power(mat, power)))
        powers1.append(var_dist(np.linalg.matrix_power(mat1, power)))
    ax1.scatter(range(1, max_pow+1), powers)
    ax1.scatter(range(1, max_pow+1),powers1)
    plt.show()
    return powers, powers1


def __main__():
    #Excedances:
    mat = DFT_excedances(4)[2]
    mat = normalize(mat, total_exced(4))
    plot_var_dist(mat, 20)
    #print(DFT_excedances(6)) -- nonzero matrix index is 4

    #Major Index:
    mat = DFT_major_index(5)[2]
    mat = normalize(mat, total_exced(5))
    mat1 = DFT_major_index(5)[4]
    mat1 = normalize(mat1, total_exced(5))
    plot_var_dist2(mat, mat1, 13)

    #Length:
    mat = DFT_length(6)[4]
    mat = normalize(mat, total_length(6))
    mat1 = DFT_length(6)[6]
    mat1 = normalize(mat1, total_length(6))
    plot_var_dist2(mat, mat1, 20)



__main__()

