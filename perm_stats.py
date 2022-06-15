from permutation import Permutation
from itertools import combinations 
from tableau import Tableau
from typing import Callable
from misc import perm_pow
import math

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

def descent(sigma: Permutation, n: int) -> int:
    """
    returns the number of descents of a given permutation
    sigma: permutation over n integers
    """
    count = 0
    for index in range(1, n):
        if sigma(index) > sigma(index+1):
            count += 1
    return count


def length(sigma: Permutation, n: int) -> int:
    """
    returns the length statistic (number of inversions) for a permutation
    note--create this function so we can use it in our analysis
    """
    return sigma.inversions()

def fixed_points(sigma: Permutation, n: int) -> int:
    '''
    Returns the number of fixed points in a permutation
    '''
    count = 0
    for i in range(1, n+1):
        if sigma(i) == i:
            count += 1
    return count

def power_function(function, k):
    return (lambda sigma, n: function(perm_pow(sigma, k), n))

def cycle_indicator(sigma: Permutation, n: int) -> int:
    generator = Permutation(*([i for i in range(2, n+1)] + [1]))
    current = generator
    for power in range(n):
        if sigma == current:
            return math.factorial(n-2) * (n*(n-1)/2  - power)
        current = current * generator
    return 0
        

def distance_from_standard(partition: tuple, n: int) -> int:
    ''' 
    Returns a permutation statistic that puts the permutation into a tableau and returns the distance the tableau is from being standard
    '''
    return (lambda sigma, n: Tableau.perm_to_tableau(sigma, partition).dist_from_standard())

def count_occurrences(sigma: Permutation, tau: Permutation, n: int, k: int) -> int:
    """
    Arguments: sigma, a permutation in Sn
               tau, a permutation in Sk
               k <= n
    Returns: the number of occurrences of tau in sigma
    """
    return

def inc_seq_k(sigma, k, n):
    one_n = list(range(1,n+1))
    count = 0
    for indices in combinations(one_n, k):
        for images in combinations(one_n, k):
            if w_gen(images, indices, sigma):
                count += 1
    return count

def w_ij(i: int, j: int) -> Callable[[Permutation, int], int]:
    """
    returns a w_ij function
    the returned function takes in a permutation sigma, and size n
    and then returns 1 if sigma maps j to i, and 0 if not
    """
    return (lambda sigma, n: int(sigma(j) == i))

def generate_w_ij(n: int) -> list:
    result = []
    for i in range(1, n+1):
        for j in range(1, n+1):
            result.append(w_ij(i, j))
    return result

def w_ij_kl(i: int, j: int, k: int, l: int) -> Callable[[Permutation, int], int]:
    """
    returns a w_ij_kl function
    the returned function takes in a permutation sigma, and size n
    and returns 1 if sigma maps k to i and l to j, and 0 otherwise
    """
    return (lambda sigma, n: int(sigma(k) == i and sigma(l) == j))

def w_gen(pattern: list, indices: list, sigma):
    for i in range(len(indices)):
        if sigma(indices[i]) != pattern[i]:
            return False
    return True

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

def normal(f: Callable[[Permutation, int], int], n: int) -> Callable:
    """
    f: permutation statistic
    Returns a function that, when called on a permutation returns f of that permutation divided by the sum of f over Sn
    """
    return (lambda sigma, n: f(sigma, n)/total(f, n))
    
