from permutation import Permutation
import math
import numpy as np

def wolfram_syntax(mat: np.array) -> None:
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
