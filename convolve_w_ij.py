import math

from perm_stats import w_ij
from perm_stats import w_ij_kl
from permutation import Permutation

# The functions in this file were written to test the hypothesis of convolving wij functions

def generate_w_ij_elem(i: int, j: int, n: int) -> dict:
    """
    Generates the element w_ij as an element in the group algebra. 
    The dictionary returned represents the group algebra element. The
    keys are permutations and the values are their coefficients.
    """
    w_ij_func = w_ij(i,j)
    w_ij_dict = {}
    scale = 1/math.factorial(n-1)
    for sigma in Permutation.group(n):
        w_ij_dict[sigma] = w_ij_func(sigma, n) * scale
    return w_ij_dict

def generate_w_ij_kl_elem(i: int, j: int, k:int, l:int, n: int) -> dict:
    """
    Generates the element w_ij_kl as an element in the group algebra. 
    The dictionary returned represents the group algebra element. The
    keys are permutations and the values are their coefficients.
    """
    w_ij_kl_func = w_ij_kl(i,j,k,l)
    w_ij_kl_dict = {}
    sn = Permutation.group(n)
    scale = 1/math.factorial(n-2)
    for sigma in sn:
        w_ij_kl_dict[sigma] = w_ij_kl_func(sigma, n) * scale
    return w_ij_kl_dict

def sum_w_ij_kl(i: int, l: int, n: int) -> dict:
    """
    Returns the sum of w_i alpha_beta l where alpha is 
    anything but i and beta is anything but l
    """
    sums = {}
    scale = 1/((n-1)**2)
    for sigma in Permutation.group(n):
        sums[sigma] = 0
    for alpha in range(1,n+1): 
        if alpha != i:
            for beta in range(1, n+1):
                if beta != l:
                    w_i_alpha_l_beta = generate_w_ij_kl_elem(i, alpha, beta, l, n)
                    for sigma in Permutation.group(n):
                        sums[sigma] += w_i_alpha_l_beta[sigma] * scale
    return sums

def convolve(i: int, j: int, k: int, l: int, n: int) -> dict:
    """
    Takes in parameters to create group elements w_ij, w_kl and
    convolves these group elements and returns a dict with the 
    resulting coefficients
    """
    wij = generate_w_ij_elem(i,j,n)
    wkl = generate_w_ij_elem(k,l,n)
    product = {}
    count = 0
    for sigma in Permutation.group(n):
        product[sigma] = 0
    for sigma in Permutation.group(n):
        for tau in Permutation.group(n):
            pi = sigma * tau
            product[pi] = product[pi] + wij[sigma]*wkl[tau]
    return product

def compare_sum_product(i: int, j: int, k: int, l: int, n) -> bool:
    """
    function for testing hypthesis for small n
    """
    prod = convolve(i,j,k,l,n)
    if j == k:
        wijkl_sum = generate_w_ij_elem(i, l, n)
    else:
        wijkl_sum = sum_w_ij_kl(i, l, n)
    equal = True
    sn = Permutation.group(n)
    for sigma in sn:
        if abs(prod[sigma] - wijkl_sum[sigma]) > 10**(-10): # Avoid computational error
            return False
    return True
    
def print_group_alg_elem(w: dict, n: int) -> None:
    """
    print function for elements in the group algebra
    prints a plus at the end which can be ignored and ends in a line break
    """
    sn = Permutation.group(n)
    for sigma in sn:
        print(str(w[sigma]) + " " + str(sigma) + " + ", end = "")
    print("\n")
