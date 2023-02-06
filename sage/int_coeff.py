# This file is dedicated to investigating whether the 1ij functions have integer coefficients when written in terms of the longest increasing subsequence basis
from math import factorial
import numpy as np
from sage.combinat.permutation import Permutation
from sage.matrix.constructor import Matrix
from sage_perm_stats import oij
from sage_misc import function_to_vector
from sage_misc import function_to_list
from itertools import product

n = 4

def dot_product(v1: Matrix, v2: Matrix) -> str:
    r = v1.nrows()
    result = ''
    for row in range(r):
        result = result + str(v1[row]) + str(v2[row]) + '+'
    return result
# 1-local T

T = Matrix([function_to_list(oij([i], [j]), n) for (i, j) in product(range(1, n+1), range(1, n+1))])




# 1-local
m = Matrix(factorial(n), (n - 1)**2 + 1 )
sigma = Permutation(list(range(1, n+1)))
col = 0
index = []
for row in range(m.nrows()):
    if sigma.longest_increasing_subsequence_length() >= n - 1:
        print(sigma, oij([1], [2])(sigma, n))
        index.append(sigma)
        m[row, col] = 1
        col += 1
    sigma = sigma.next()
# print(m)

v_index = Matrix(index)
m_proj = T * m


for i in range(1, n+1):
    for j in range(1, n+1):
        b = function_to_vector(oij([i], [j]), n)
        b_proj = T * b
        print(i, ' -> ', j)
        lc = m_proj.solve_right(b_proj)
        print(dot_product(lc, v_index))

        

# 2-local
