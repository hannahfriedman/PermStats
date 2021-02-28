import numpy as np
from numpy.linalg import matrix_power
import math
import copy
from Tableau import Tableau
from permutation import Permutation
#Docs for Permutation class: 
#https://permutation.readthedocs.io/en/stable/_modules/permutation.html#Permutation.cycle

def fac(n):
    if n == 0:
        return 1
    else:
        return n*fac(n-1)

def DTF(n, f):
    """
    n--int n in S_n
    f--dict function from S_n to Z
    """
    return

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
    nfac = fac(n)
    
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
    rev_partitions = partitions.reverse()
    tableaux_by_shape = [tableaux_shape(n, partition) for partition in partitions]
    rho = {}
    for i in range(1,n):
        representation = []
        for shape in tableaux_by_shape:
            sort_tableaux(n, shape)
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


def generate_tableaux(n):
    '''
    Generates all tableaux of a given size n
    '''
    if n == 1:
        return [Tableau([[1]])]
    else:
        prev_tab = generate_tableaux(n-1)
        ans = []
        for tab in prev_tab:
            for i in range(len(tab.data)):
        # Aldrin: this first if statement doesn't seem necessary so it is removed for now
                # if tab.size == 1: 
                    # ans += [Tableau([tab.data[i] + [n]])]
        # adds n to end of row if allowed
                if i == 0 or len(tab.data[i]) < len(tab.data[i-1]):
                    ans += [Tableau(tab.data[0:i] + [tab.data[i] + [n]] + tab.data[i+1:])]
        # adds n as a new row
        ans += [Tableau(tab.data + [[n]]) for tab in prev_tab]

        return ans



def tableaux_shape(n, partition):
    '''
    Generates all tableaux of size n, that fit a given partition
    '''
    ans = []
    tableaux_list = generate_tableaux(n)
    for tab in tableaux_list:
        if tab.shape == partition:
            ans += [tab]
    return ans
      
def sort_tableaux(n, tableaux):
    switched = True
    while switched:
        switched = False
        for i in range(len(tableaux) - 1):
            if tableaux[i] < tableaux[i+1]:
                tableaux[i], tableaux[i+1] = tableaux[i+1], tableaux[i]
                switched = True
    return

