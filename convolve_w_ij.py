from perm_stats import w_ij
from perm_stats import w_ij_kl
from permutation import Permutation
import math

def generate_w_ij_elem(i: int, j: int, n: int) -> dict:
    """
    Generates the element w_ij as an element in the group algebra. 
    The dictionary returned represents the group algebra element. The
    keys are permutations and the values are their coefficients.
    """
    w_ij_func = w_ij(i,j)
    w_ij_dict = {}
    sn = Permutation.group(n)
    scale = 1/math.factorial(n-1)
    for sigma in sn:
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
    Returns the sum of w_i alphas and w_beta ls where alpha is 
    anything by i and beta is anything by l
    """
    sums = {}
    sn = Permutation.group(n)
    scale = 1/(n-1)**2
    for sigma in sn:
        sn = Permutation.group(n)
        sums[sigma] = 0
    for alpha in range(1,n+1): 
        if alpha != i:
            for beta in range(1, n+1):
                if beta != l:
                    w_i_alpha_l_beta = generate_w_ij_kl_elem(i, alpha, beta, l, n)
                    sn1 = Permutation.group(n)
                    for sigma in sn1:
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
    sn0 = Permutation.group(n)
    sn1 = Permutation.group(n)
    count = 0
    for sigma in sn0:
        product[sigma] = 0
    for sigma1 in sn1:
        sn2 = Permutation.group(n)
        for sigma2 in sn2:
            product[sigma1 * sigma2] = product[sigma1 * sigma2] + wij[sigma1]*wkl[sigma2]
            if wij[sigma1]*wkl[sigma2] != 0:
                print(str(sigma1 * sigma2) + " ", end = '')
            # if sigma1 * sigma2 == Permutation(2,3,1,4):
            #     count += 1
            #     print(str(wij[sigma1]*wkl[sigma2]) + "*" + str(sigma1) + "*" + str(sigma2) + ": " + str(count))
        print("")
    return product
    
    
def print_group_alg_elem(w: dict, n: int) -> None:
    """
    print function for elements in the group algebra
    prints a plus at the end which can be ignored and ends in a line break
    """
    sn = Permutation.group(n)
    for sigma in sn:
        print(str(w[sigma]) + " " + str(sigma) + " + ", end = "")
    print("\n")



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

# Test hypothesis for small n; do not try to test for n>5 
# n = 5
# for i in range(1,n+1):
#     for j in range(1,n+1):
#         for k in range(1,n+1):
#             for l in range(1,n+1):
#                 if (not compare_sum_product(i,j,k,l,n)):
#                     print("whoof")


