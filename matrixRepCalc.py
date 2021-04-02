import numpy as np
from numpy.linalg import matrix_power
import math
import copy
from Tableau import Tableau
from excedances import count_excedances
from permutation import Permutation
from majorIndex import calc_major_index
from misc import matlab_syntax
from wij import w_ij

#Docs for Permutation class: 
#https://permutation.readthedocs.io/en/stable/_modules/permutation.html#Permutation.cycle

def DFT_length(n):
    """
    n--int n, size of permutation group S_n
    f--dict function from S_n to Z
    """
    nfac = math.factorial(n)
    sn = Permutation.group(n)
    dft = []
    rho = matrix_rep(n)
    for i in range(nfac):
        perm = next(sn)
        length = perm.inversions()
        if i == 0:
            for mat in rho[perm]:
                dft.append(length*mat)
        else:
            for i in range(len(rho[perm])):
                dft[i] = np.add(dft[i], length*rho[perm][i])
    return adjust_zeros(dft)

def DFT_excedances(n):
    """
    n--int n in S_n
    """
    nfac = math.factorial(n)
    sn = Permutation.group(n)
    dft = []
    rho = matrix_rep(n)
    for i in range(nfac):
        perm = next(sn)
        excedances = count_excedances(perm, n)
        if i == 0:
            for mat in rho[perm]:
                dft.append(excedances*mat)
        else:
            for i in range(len(rho[perm])):
                dft[i] = np.add(dft[i], excedances*rho[perm][i])

    return adjust_zeros(dft)


def DFT_major_index(n):
    """
    n--int n in S_n
    """
    nfac = math.factorial(n)
    sn = Permutation.group(n)
    dft = []
    rho = matrix_rep(n)
    for i in range(nfac):
        perm = next(sn)
        major_index = calc_major_index(perm, n)
        if i == 0:
            for mat in rho[perm]:
                dft.append(major_index*mat)
        else:
            for i in range(len(rho[perm])):
                dft[i] = np.add(dft[i], major_index*rho[perm][i])

    return adjust_zeros(dft)

def DFT_w_ij(n, i, j):
    """
    n--int n in S_n
    """
    nfac = math.factorial(n)
    sn = Permutation.group(n)
    dft = []
    rho = matrix_rep(n)
    for index in range(nfac):
        perm = next(sn)
        wij = w_ij(perm, j, i)
        if index == 0:
            for mat in rho[perm]:
                dft.append(wij*mat)
        else:
            for rep_index in range(len(rho[perm])):
                dft[rep_index] = np.add(dft[rep_index], wij*rho[perm][rep_index])
    return adjust_zeros(dft)

def adjust_zeros(dft):
    for mat in dft:
        for row in range(mat.shape[0]):
            for col in range(mat.shape[1]):
                if abs(mat[row, col]) <= 10**-10:
                    mat[row,col] = 0
    return dft

def matrix_rep(n):
    """
    n--int n in S_n
    returns a dict that maps the 2-cycles of S_n to their orthogonal matrix representations
    """
    rho_gen = matrix_rep_gen(n)
    rho = {}
    sn = factor_sn(n)
    for perm in sn:
        key = Permutation.cycle()
        val = []
        for mat in rho_gen[Permutation.cycle(1,2)]:
            val.append(matrix_power(mat, 2))
        for transposition in perm:
            if transposition != Permutation.cycle():
                key *= transposition
                for i in range(len(val)):
                    val[i] = np.matmul(val[i], rho_gen[transposition][i])
        rho[key] = val
    return rho


def factor_sn(n):
    """
    n--int n in S_n
    returns a list of lists of factorizations of all elements of Sn
    """
    if n == 2:
        return [[Permutation.cycle()], [Permutation.cycle(1,2)]]
    else:
        sn_minus_one = factor_sn(n-1)
        sn = sn_minus_one
        prev_coset = sn_minus_one
        curr_coset = []
        for i in range(1, n):
            for perm in prev_coset:
                newPerm = perm + [Permutation.cycle(n-i,n-i+1)]
                curr_coset.append(newPerm)
            sn = sn + curr_coset
            prev_coset = curr_coset
            curr_coset = []
        return sn



def matrix_rep_transpositions(n):
    """
    n--int n in S_n
    returns a dict that maps the 2-cycles of S_n to their orthogonal matrix representations
    With help from formulas obtained in: https://math.stackexchange.com/questions/3420570/writing-permutations-as-products-of-adjacent-transposition 
    """
    rho_gen = matrix_rep_gen(n)
    rho = rho_gen
    
    # Calcultes the matrix representation for all 2-cycles
    for diff in range(2, n):
        for startVal in range(1, n - diff + 1):
            endVal = startVal + diff
            permutation = Permutation.cycle(startVal, endVal)
            perm_factor1 = Permutation.cycle(endVal - 1, endVal)
            perm_factor2 = Permutation.cycle(startVal, endVal-1)
            matrix_factor1 = rho[perm_factor1]
            matrix_factor2 = rho[perm_factor2]
            matrix_rep = [np.matmul(np.matmul(matrix_factor1[i], 
                                    matrix_factor2[i]), 
                                    matrix_factor1[i])
                                    for i in range(len(matrix_factor1))]
            rho[permutation] = matrix_rep
    
    # Matrix representation for the rest
    return rho

        

def matrix_rep_gen(n):
    """
    n--int n in S_n
    returns a dict that maps the generators of S_n to their orthogonal matrix representations
    """
    partitions = generate_partitions(n)
    partitions.reverse()
    #print(partitions)
    tableaux_by_shape = [tableaux_shape(n, partition) for partition in partitions]
    rho = {}
    for i in range(1,n):
        representation = []
        for shape in tableaux_by_shape:
            sort_tableaux(n, shape)
            print(shape)
            rep = np.zeros((len(shape), len(shape)))
            for index in range(len(shape)):
                tableau = Tableau(shape[index].data)
                rep[index, index] = 1/(tableau.signed_distance(i))
                switched = tableau.switch(i)
                if switched.is_standard():
                    switched_index = 0
                    for j in range(len(shape)):
                        if shape[j] == switched:
                            switched_index = j
                            break
                    rep[switched_index, index] = math.sqrt(1 - (shape[index].signed_distance(i))**(-2))
            representation.append(rep)
        rho[Permutation.cycle(i, i+1)] = representation
    return rho


def test_matrices(n):
    rho = matrix_rep(n)
    for i in range(1, n):
        for mat in rho["(" + str(i) + "," + str(i+1) + ")"]:
            print(matrix_power(mat, 2))
    print(40*"-")
    if n >= 4:
        for i in range(1, n-2):
            for j in range(len(rho["(" + str(i) + "," + str(i+1) + ")"])):
                mat1 = rho["(" + str(i) + "," + str(i+1) + ")"][j]
                mat2 = rho["(" + str(i+2) + "," + str(i+3) + ")"][j]
                #print(np.matmul(mat1, mat2))
                #print(np.matmul(mat2, mat1))
                print(np.matmul(mat1, mat2) == np.matmul(mat2, mat1))
    print(40*"-")
    if n > 2:
        for i in range(1, n-1):
            for j in range(len(rho["(" + str(i) + "," + str(i+1) + ")"])):
                mat1 = rho["(" + str(i) + "," + str(i+1) + ")"][j]
                mat2 = rho["(" + str(i+1) + "," + str(i+2) + ")"][j]
                LHS = np.matmul(np.matmul(mat1, mat2), mat1)
                RHS = np.matmul(np.matmul(mat2, mat1), mat2)
                print(LHS)
                print(RHS)
                #print(np.matmul(mat1, mat2) == np.matmul(mat2, mat1))
    return 

def generate_partitions(n):
    '''
    Generates all partitions of size n
    Returns a list of lists, showing the size of each partition
    '''
    ans = []
    if n == 1:
        return [[1]]
    elif n == 0:
        return [[]]
    for x in range(1, n):
        ans += [[x] + part for part in generate_partitions(n-x)]
    return remove_dubs(ans) + [[n]]

def remove_dubs(partition):
    ''' 
    Removes duplicates in a list of lists
    Makes sure that any inner lists are sorted first (treats inner lists as multi-sets)
    ''' 
    for part in partition:
        part.sort()
        part.reverse()
    result = []
    for part in partition:
        if part not in result:
            result.append(part)
    return result

def __main__():
    """
    rho = matrix_rep(3)
    print(rho)
    sn = Permutation.group(3)
    for sigma in sn:
    print(rho[sigma][1])"""
    for arr in DFT_w_ij(5,5,5):
        print(arr)

__main__()