from permutation import Permutation
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
        for col in range(len(mat)):
            string += " "
            string += str(mat[row, col])
        if row == len(mat) - 1:
            string += "]"
        else:
            string += ";"
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

def perm_pow(sigma, k):
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

def function_to_vector(function, n):
    result = []
    for g in Permutation.group(n):
        result.append(function(g, n))
    return result

def vector_to_function(vector, n):
    d = {}
    index = 0
    for g in Permutation.group(n):
        d[g] = vector[index]
        index += 1
    return (lambda sigma, n: d[sigma])

def left_action_on_vector(vector, current_basis, sigma):
    return [vector[current_basis.index(sigma * permutation)] for permutation in current_basis]

def left_right_action_on_vector(vector, current_basis, sigma, tau):
    return [vector[current_basis.index(sigma * permutation * tau)] for permutation in current_basis]

def add_functions(f1, f2, n):
    v1 = function_to_vector(f1, n)
    v2 = function_to_vector(f2, n)
    v = np.array(v1) + np.array(v2)
    return vector_to_function(v, n)

def add_weighted_functions(f1, f2, w1, w2, n):
    v1 = function_to_vector(f1, n)
    v2 = function_to_vector(f2, n)
    v = w1 * np.array(v1) + w2* np.array(v2)
    return vector_to_function(v, n)

def distribution(v):
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
        
def mats_into_vec(dim, mats):
    v = np.zeros((dim,))
    count = 0
    for m in mats:
        for row in range(m.shape[0]):
            for col in range(m.shape[1]):
                v[count] = m[row, col]
                count += 1
    return v
