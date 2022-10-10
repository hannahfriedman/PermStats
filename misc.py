from permutation import Permutation
from typing import Callable
import math
import numpy as np
from sympy import Matrix

def wolfram_syntax(mat: np.array) -> None:
    '''
    Returns a numpy arrray  as a string in wolfram alpha syntax
    '''
    string = "{"
    for row in range(len(mat)):
        string+="{"
        for col in range(len(mat[0])-1):
            string += str(mat[row, col])
            string += ","
        string += str(mat[row, len(mat[0])-1])
        if row == len(mat) - 1:
            string += "}"
        else:
            string += "},"
    string+= "}"
    print(string)

def matlab_syntax(mat: np.array) -> None:
    '''
    Returns a numpy  array  as a string in matlab syntax
    '''
    string = "["
    for row in range(len(mat)):
        for col in range(mat.shape[1]):
            string += " "
            string += str(mat[row, col])
        if row == len(mat) - 1:
            string += "]"
        else:
            string += ";\n"
    print(string)

def latex_syntax(mat: np.array) -> None:
    """
    returns a string that prints the matrix to insert into latex
    for example:
        a & b & c \\
        d & e & f \\
        g & h & i
    """
    string = ""
    for row in range(mat.shape[0]):
        for col in range(mat.shape[1]):
            string += str(mat[row, col])
            string += " & "
        string += r'\\ '   # to avoid treating \\ as special character
        string += "\n"
    string = string[:-4]   # to remove last end line
    print(string)

def adjust_zeros(dft: dict) -> dict:
    """
    Accounts for computational errors by setting very small values equal to 0
    """
    for mat in dft:
        for row in range(mat.shape[0]):
            for col in range(mat.shape[1]):
                if abs(mat[row, col]) <= 10**-10:
                    mat[row,col] = 0
    return dft

def falling_factorial(n: int, k:int) -> int:
    """
    returns n falling factorial k
    """
    product = 1
    while n > k:
        product *= n
        n -= 1
    return product

def choose(n: int, k: int) -> int:
    """
    Returns n choose k
    """
    if k > n:
        return 0
    return falling_factorial(n, k)//math.factorial(n - k)

def eig_multiplicity(m: np.matrix) -> dict:
    """
    Returns a dictionary, the keys of which are eigenvalues and the values of which are the eigenvalue's
    multiplicity
    """
    dictOfEigs = {}
    eigs = np.linalg.eig(m)[0]
    # Add each eigenvalue to the dictionary
    for eig in eigs:
        roundedEig = round(eig)
        if roundedEig not in dictOfEigs.keys():
            dictOfEigs[roundedEig] = 1
        else:
            dictOfEigs[roundedEig] += 1
    return dictOfEigs

def perm_pow(sigma: Permutation, k: int) -> Permutation:
    ''' Returns sigma^k '''
    if k == 0:
        return Permutation()
    elif k%2 == 1:
        return sigma * perm_pow(sigma, k - 1)
    else:
        half = perm_pow(sigma, k//2)
        return half *  half

def nullspace(m: np.matrix) -> list:
    spanning_set = []
    U, S, Vt = np.linalg.svd(m)
    for i in range(len(S)):
        if np.linalg.norm(S[i]) <= 10**-10:
            spanning_set.append(Vt[i]/Vt[i, 0])
    return spanning_set

def function_to_vector(function: Callable[[Permutation, int], int], n: int) -> list:
    ''' Given a function from S_n -> Z, return a list of the values of the function on S_n'''
    result = []
    for g in Permutation.group(n):
        result.append(function(g, n))
    return result

def vector_to_function(vector, n: int) -> Callable[[Permutation, int], int]:
    ''' Vector can be any iterable/vector/list that has length n!
    Returns the function version of vector
    '''
    d = {}
    index = 0
    for g in Permutation.group(n):
        d[g] = vector[index]
        index += 1
    return (lambda sigma, n: d[sigma])

def dft_to_vector(mats: list) -> np.array:
    '''
    Given a list of matrices, return all entries in a vector, row by row in the order of the original list
    '''
    len_v = sum([mat.shape[0]*mat.shape[1] for mat in mats])
    v = np.zeros((len_v))
    index = 0
    for mat in mats:
        for row in mat:
            for entry in row:
                v[index] = entry
                index += 1
    return v

def left_action_on_vector(vector: list, current_basis: list, sigma: Permutation) -> list:
    '''
    vector: a list of ints
    current_basis: a list of permutations
    together, vector and current_basis correspond to a function on S_n, given by the sum of vector[i]*current_basis[i] over i
    When we apply  sigma on the left, we get vector[i]*sigma*current_basis[i]. 
    We return the vector that puts the coefficients back in the order of current basis, i.e. v such that vector[i]*sigma*current_basis[i] = v[i] * current_basis[i]
    '''
    return [vector[current_basis.index(sigma.inverse() * permutation)] for permutation in current_basis]

def left_right_action_on_vector(vector: list, current_basis: list, sigma: Permutation, tau: Permutation) -> list:
    ''' Does the same thing as the previous function but on the left and right '''
    return [vector[current_basis.index(sigma.inverse() * permutation * tau.inverse())] for permutation in current_basis]

def distribution(v) -> list:
    '''
    v: vector/list/tuple
    returns a vector p such that p[i] is the number of times i appears in v
    '''
    result = []
    for i in v:
        while i >= len(result):
            result.append(0)
        result[i] += 1
    return result

def similar(m1, m2):
    ''' Return P such that m2 = Pm1P^(-1) or false if m1, m2 are not similar'''
    m1 = Matrix(m1)
    m2 = Matrix(m2)
    p1, j1 = m1.jordan_form()
    p2, j2 = m2.jordan_form()
    p1 = np.array(p1).astype(np.float64)
    j1 = np.array(j1).astype(np.float64)
    p2 = np.array(p2).astype(np.float64)
    j2 = np.array(j2).astype(np.float64)
    if not (j1 == j2).all():
        return False
    return p2 @ np.linalg.inv(p1)

def function_sum_to_vector(n, *args):
    result = np.zeros((1, math.factorial(n)))
    for f in args:
        result += np.array(function_to_vector(f, n))
    return result

def function_sum(a, b, a_coeff=1, b_coeff=1):
    ''' Adds two permutation statistics with coefficients'''
    if a == None:
        return (lambda sigma, n: b_coeff * b(sigma, n))
    if b == None:
        return (lambda sigma, n: a_coeff * a(sigma, n))
    return (lambda sigma, n: a_coeff * a(sigma, n) + b_coeff * b(sigma, n))

def mats_into_vec(dim, mats):
    '''
    Reshapes a list of matrices into a vector
    '''
    v = np.zeros((dim,))
    count = 0
    for m in mats:
        for row in range(m.shape[0]):
            for col in range(m.shape[1]):
                v[count] = m[row, col]
                count += 1
    return v

def put_in_relative_order(nums, order):
    '''
    Takes in two lists of equal length and returns the contents of numbers in the relative order specified by order
    '''
    if len(nums) != len(order):
        raise Exception("Inputs must have the same length.")
    nums = list(nums)
    order = list(order)
    result = [0] * len(nums)
    while nums != []:
        i = order.index(max(order))
        result[i] = max(nums)
        nums.pop(nums.index(max(nums)))
        order[i] = -1
    return result
        
        
def apply_permutation(values: list, sigma: Permutation):
    '''
    values is an iterable with lists of integers, i.e. [(1, 2, 3), (4, 6), (5,)]
    '''
    return [[sigma(k) for k in val] for val in values]
    

def less_than_set(set1: set, set2: set) -> bool:
    ''' returns true if set1 < set2
    set1 < set2 if min(set1) < min(set2). If they are equal, use recursion.'''
    if len(set1) != len(set2):
        raise Exception("Sets have different sizes.")
    l1 = sorted(list(set1))
    l2 = sorted(list(set2))
    for i in range(len(l1)):
        if l1[i] >= l2[i]:
            return False
    return True
    # if set1 == set2:
    #     return False
    # if min(set1) == min(set2):
    #     new1 = set1.copy()
    #     new1.remove(min(set1))
    #     new2 = set2.copy()
    #     new2.remove(min(set2))        
    #     return less_than_set(new1, new2)
    # return sum(set1) < sum(set2)
    # return min(set1) < min(set2)
