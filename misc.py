from permutation import Permutation
import math
def wolfram_syntax(mat):
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

def matlab_syntax(mat):
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

def latex_syntax(mat):
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

def normalize(mat, stat):
    for row in range(mat.shape[0]):
        for col in range(mat.shape[1]):
            mat[row, col] =  mat[row, col]/stat
    return mat


def adjust_zeros(dft):
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