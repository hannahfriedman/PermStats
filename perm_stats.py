from permutation import Permutation
from tabloid import Tabloid
from itertools import combinations 
from tableau import Tableau
from typing import Callable
from misc import perm_pow
from misc import put_in_relative_order
from misc import function_sum
import math

################
## Statistics ##
################

def excedances(sigma: Permutation, n: int) -> int:
    """
    returns the number of excedances of a given permutation
    sigma: given permutation over n integexrs
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

def total_displacement(sigma: Permutation, n: int) -> int:
    '''
    returns the total displacement of a permutation, i.e. |sigma(i) - i| summeder over i = 1, ..., n+1
    sigma: permutation of n integers
    '''
    return sum([abs(sigma(i) - i) for i in range(1, n+1)])


def length(sigma: Permutation, n: int) -> int:
    """
    returns the length statistic (number of inversions) for a permutation
    """
    return sigma.inversions()

def fixed_points(sigma: Permutation, n: int) -> int:
    '''
    Returns the number of fixed points in a permutation
    sigma: permutation over n integers
    '''
    count = 0
    for i in range(1, n+1):
        if sigma(i) == i:
            count += 1
    return count


def pairs_of_fixed_points(sigma: Permutation, n: int) -> int:
    ''' 
    Returns the number of pairs of fixed points in a permutation, i.e. the number of ordered pairs (i,j) s.t. sigma(i) = i and sigma(j) = j
    sigma: permutation over n integers
    '''
    count = 0
    for i, j in combinations(list(range(1, n+1)), 2):
        count += 2 * int(sigma(i) == i and sigma(j) == j)  # Multiply by 2 to count (i,j) and (j,i)
    return count

def unordered_fixed_points(sigma: Permutation, n: int) -> int:
    ''' 
    Returns the number of 2-sets that sigma maps to itself, i.e. the number 2-sets {i,j} with either
        sigma(i) = i, sigma(j) = j
    or  sigma(i) = j, sigma(j) = i
    sigma: permutation over n integers
    '''
    count = 0
    for i, j in combinations(list(range(1, n+1)), 2):
        count += int(sigma(i) == i and sigma(j) == j)
        count += int(sigma(i) == j and sigma(j) == i)        
    return count

def position_of_n(sigma: Permutation, n:int) -> int:
    '''
    Returns i such that sigma(i) = n
    sigma: Permutation on n integers
    '''
    return sigma.inverse()(n)

def cycle_indicator(sigma: Permutation, n: int) -> int:
    '''
    Returns (n-2)!(n*(n-1)/2 - k - 1) if sigma = (2 3 ... n 1)^k and zero otherwise
    This function has the same one-local projection as excedances.
    '''
    generator = Permutation(*([i for i in range(2, n+1)] + [1]))
    current = generator
    for power in range(n):
        if sigma == current:
            return math.factorial(n-2) * (n*(n-1)/2  - power)
        current = current * generator
    return 0

def all_ones(sigma: Permutation, n: int) -> int:
    return 1

########################################
## Functions that Generate Statistics ##
########################################

def count_cycles(k: int) -> Callable[[Permutation, int], int]:
    ''' 
    Returns a function that counts the number of k-cycles in a permutation
    sigma: Permutation on n integers
    '''
    def f(sigma: Permutation, n: int) -> int:
        '''
        Return the number of k-cycles in sigma, a permutation on n elements
        '''
        count = 0
        for sets in combinations(list(range(1, n+1)), k):
            current = sets[0]                  # Pick an arbitrary  element in the k-cycle to start the cycle
            cycle_length = 0                   # Keep track of the cycle length
            while current in sets:             # If we map outside of the current k-set, go to the next k-set
                current = sigma(current)       # Update the value and cycle length
                cycle_length += 1
                if current == sets[0]:         # If we have completed the cycle and it is the correct length,
                    if cycle_length == k:      # count the cycle, otherwise move on to the next cycle
                        count += 1
                    else:
                        break
        return count
    return f

def power_function(function: Callable[[Permutation, int], int], k: int) -> Callable[[Permutation, int], int]:
    '''
    Given a function f, return a new function that returns f(sigma^k)
    '''
    return (lambda sigma, n: function(perm_pow(sigma, k), n))


def distance_from_standard(partition: tuple, n: int, method=0) -> Callable[[Permutation, int], int]:
    ''' 
    Returns a permutation statistic that puts the permutation into a tableau left to right and top to bottom and returns the distance the tableau is from being standard
    '''
    return (lambda sigma, n: Tableau.perm_to_tableau(sigma, partition).dist_from_standard(method=method))

def count_occurrences(tau: Permutation, n: int, k: int) -> Callable[[Permutation, int], int]:
    """
    Arguments: sigma, a permutation in Sn
               tau, a permutation in Sk
               k <= n
    Returns: a function counting the number of occurrences of tau in sigma
    """
    one_n = list(range(1,n+1))
    order = [tau(i) for i in range(1, k+1)]
    def f(sigma, n):
        result = 0
        for indices in combinations(one_n, k):
            for images in combinations(one_n, k):
                image = put_in_relative_order(images, order)
                result += pattern_indicator(sigma, indices, image)
        return result
    return f

def inc_seq_k(k: int) -> int:
    ''' Return the number occurances of an increasing sequences of length k in sigma '''
    def f(sigma: Permutation, n: int) -> int:
        one_n = list(range(1,n+1))
        count = 0
        for indices in combinations(one_n, k):
            for images in combinations(one_n, k):
                if pattern_indicator(sigma, indices, images):
                    count += 1
        return count
    return f

def share_fixed(tab: Tabloid, n: int):
    ''' Given a tabloid, return the associated function '''
    def f(t: Tabloid, n: int):
        # count = 0
        # for row in range(len(tab.data)):
        #     # count += len(tab.data[row].intersection(t.data[row]))
        #     # count += len(tab.data[row]) * int(tab.data[row] == t.data[row])
        #     # count += int(tab.data[row].intersection(t.data[row]) != set())
        # return count
        return int(t == tab)
    return f

def tabloid_to_perm_stat(tab: Tabloid):
    ''' Given a tabloid, use cosets of its stabilizer to pull the share_fixed function back to the group'''
    def g(sigma: Permutation, n: int):
        return share_fixed(tab, n)(tab.apply_permutation(sigma), n)
    return g

def average_tabloid_statistic(shape: tuple, n: int):
    ''' Average functions on the group generated by using stabilizers of all elements of the set'''
    stat = None
    for tab in Tabloid.generate_by_shape(shape, n):
        stat = function_sum(stat, tabloid_to_perm_stat(tab))
    return function_sum(None, stat)#, b_coeff=1/math.factorial(n))

#########################
## Indicator Functions ##
#########################

def w_ij(i: int, j: int) -> Callable[[Permutation, int], int]:
    """
    returns a w_ij function
    the returned function takes in a permutation sigma, and size n
    and then returns 1 if sigma maps j to i, and 0 if not
    """
    return (lambda sigma, n: int(sigma(j) == i))

def generate_w_ij(n: int) -> list:
    ''' Generate a list of all the wij functions for particular n in the order w11, w12, ..., w1n, w21, ..., wnn '''
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

def w_ij_kl_unordered(i: int, j: int, k: int, l:int) -> Callable[[Permutation, int], int]:
    ''' Return true if {k,l} -> {i,j} '''
    return (lambda sigma, n: int((sigma(k) == i and sigma(l) == j) or (sigma(l) == i and sigma(k) == j)))

def pattern_indicator(sigma: Permutation, indices: list, images: list) -> bool:
    ''' Return true if sigma(indices) = images in order, false otherwise'''
    for i in range(len(indices)):
        if sigma(indices[i]) != images[i]:
            return False
    return True

def w_I_J(I: list, J: list) -> Callable[[Permutation, int], int]:
    '''
    I, J can be lists of tuples or lists of integers
    Additionally, the length of I[i] must equal the length of J[i] for all i.
    Returns a function that returns 1 if, for all i, (J[i]) = I[i] and 0 otherwise
    '''
    for i in range(len(I)):
        if len(I[i]) != len(J[i]):
            raise Exception("Shapes of I and J don't match")
    def f(sigma: Permutation, n: int) -> int:
        count = 0
        for i in range(len(I)):
            for j in J[i]:
                if sigma(j) not in I[i]:
                    return 0
        return 1
    return f

###########################
## Normalizing Functions ## 
###########################

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
    

