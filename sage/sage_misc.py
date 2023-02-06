from sage.combinat.permutation import Permutation
from sage.matrix.constructor import Matrix
from math import factorial

def function_to_vector(func, n):
    sigma = Permutation(list(range(1, n+1)))
    v = Matrix(factorial(n), 1)
    for row in range(v.nrows()):
        v[row] = func(sigma, n)
        sigma = sigma.next()
    return v


def function_to_list(func, n):
    sigma = Permutation(list(range(1, n+1)))
    v = []
    for _ in range(factorial(n)):
        v.append(func(sigma, n))
        sigma = sigma.next()
    return v
