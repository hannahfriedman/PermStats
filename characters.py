import numpy as np
from permutation import Permutation
from dft import dft
from dft import dft_natural
from dft import one_local_dft_natural
from misc import function_to_vector
from misc import add_weighted_functions
from perm_stats import fixed_points
from perm_stats import w_ij
from perm_stats import excedances
from perm_stats import length
from math import factorial
def chi_one(sigma, n):
    return 1

def chi_two(sigma, n):
    cycles = sigma.to_cycles()
    if len(cycles) == 0:
        return 3
    elif len(cycles) == 2:
            return -1
    elif len(cycles) == 1:
        cycle = cycles[0]
        if len(cycle) == 2:
            return 1
        elif len(cycle) == 3:
            return 0
        else:
            return -1


n = 171
print(one_local_dft_natural([(1, 1, 1), (2, 2, 1), (3, 3, 1), (4, 4, 1), (5, 5, 1)], n))
# print(function_to_vector(chi_two, n))
# print(dft(chi_one, n))
# print(dft_natural(add_weighted_functions(chi_one, chi_two, 1, 1, n), n))
# print(function_to_vector(add_weighted_functions(chi_one, chi_two, 1, 1, n), n))
# print(function_to_vector(fixed_points, n))
# print(dft_natural(w_ij(1, 3), n))
# for n in range(5, 6):
#     m = dft_natural(length, n)
#     print(m)
    # print(np.round_(np.linalg.eig(m)[0], 5))

