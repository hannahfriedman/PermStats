from permutation import Permutation
from typing import Any

def excedances(sigma: Permutation, n: int) -> int:
    """
    returns the number of excedances of a given permutation
    sigma: given permutation over n integers
    """
    count = 0
    for i in range(1, n+1):
        if i < sigma(i):
            count += 1
    return count

def major_index(sigma: Permutation, n: int) -> int:
    """
    returns the major index of a given permutation
    sigma: permutation over n integers
    """
    count = 0
    for index in range(1, n):
        if sigma(index) > sigma(index+1):
            count += index
    return count

def length(sigma: Permutation, n: int) -> int:
    """
    returns the length statistic (number of inversions) for a permutation
    note--create this function so we can use it in our analysis
    """
    return sigma.inversions()

def w_ij(i: int, j: int) -> Any:
    """
    returns a w_ij function
    the returned function takes in a permutation sigma, and size n
    and then returns 1 if sigma maps j to i, and 0 if not
    """
    return (lambda sigma, n: int(sigma(j) == i))

def w_ij_kl(i: int, j: int, k: int, l: int) -> Any:
    """
    returns a w_ij_kl function
    the returned function takes in a permutation sigma, and size n
    and returns 1 if sigma maps k to i and l to j, and 0 otherwise
    """
    return (lambda sigma, n: int(sigma(k) == i and sigma(l) == j))