from permutation import Permutation
from typing import Callable

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

def count_occurrences(sigma: Permutation, tau: Permutation, n: int, k: int) -> int:
    """
    Arguments: sigma, a permutation in Sn
               tau, a permutation in Sk
               k <= n
    Returns: the number of occurrences of tau in sigma
    """
    return



def w_ij(i: int, j: int) -> Callable[[Permutation, int], int]:
    """
    returns a w_ij function
    the returned function takes in a permutation sigma, and size n
    and then returns 1 if sigma maps j to i, and 0 if not
    """
    return (lambda sigma, n: int(sigma(j) == i))

def w_ij_kl(i: int, j: int, k: int, l: int) -> Callable[[Permutation, int], int]:
    """
    returns a w_ij_kl function
    the returned function takes in a permutation sigma, and size n
    and returns 1 if sigma maps k to i and l to j, and 0 otherwise
    """
    return (lambda sigma, n: int(sigma(k) == i and sigma(l) == j))

def total(f: Callable[[Permutation, int], int], n: int) -> int:
    """
    Returns the sum of f called on every permutation in Sn
    f should take inputs sigma (permutation) and n (int) and return an int
    """
    count = 0
    sn = Permutation.group(n)
    for sigma in sn:
        count += f(sigma, n)
    return count