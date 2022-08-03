from perm_stats import w_I_J
from perm_stats import w_ij
from perm_stats import w_ij_kl
from perm_stats import descent
from perm_stats import length
from perm_stats import distance_from_standard
from perm_stats import average_tabloid_statistic
from perm_stats import share_fixed
import numpy  as np
from dft import dft
from misc import apply_permutation
from misc import function_sum
from permutation import Permutation
from math import factorial
from itertools import product
from tabloid import Tabloid
from dft import sort_partitions
from dft import generate_partitions

n = 4
# t = Tabloid([[3, 4], [1], [2]])
for tab in Tabloid.generate_by_shape((2, 1, 1), n):
    print(share_fixed(tab.apply_permutation(Permutation(1, 2, 4, 3)), n)(tab, n), share_fixed(tab.apply_permutation(Permutation(2, 1)), n)(tab, n))    

# print(t.descents())
n = 4
for sigma in Permutation.group(n):
    print(sigma, round(average_tabloid_statistic((1,1,1,1), n)(sigma, n)*24, 3), round(average_tabloid_statistic((2,1,1), n)(sigma, n)*24, 3), round(average_tabloid_statistic((3,1), n)(sigma, n)*24, 3))
n = 5
for p in sort_partitions(generate_partitions(n)):
    if len(p) == 1:
        print(p)
        print(dft(average_tabloid_statistic(p, n), n))
    else:
        do_print = True
        # for i in range(1, len(p)):
        #     if p[i] != 1:
        #         do_print = False
        if do_print:
            print(p)
            print([m[0,0]*m.shape[0]/120 for m in dft(average_tabloid_statistic(p, n), n)])

k = 2
n = 5
I_0 = [(1, 2), (3, 4, 5)]

f = None
# for sigma in Permutation.group(n):
#     f = function_sum(f, w_I_J(apply_permutation(I_0, sigma), I_0), b_coeff=(sum(I_0[1]) - sum(apply_permutation(I_0, sigma)[1]))/(factorial(len(I_0[0])) * factorial(len(I_0[1]))))

# print(f(Permutation(), n))
# print(f(Permutation(1, 2, 4, 3, 5), n))
# print(f(Permutation(1, 2, 5, 4, 3), n))
# print(f(Permutation(1, 4, 3, 2, 5), n))
# print(f(Permutation(1, 5, 3, 4, 2), n))
# print(f(Permutation(1, 4, 5, 2, 3), n))
# print(f(Permutation(4, 2, 3, 1, 5), n))
# print(f(Permutation(5, 2, 3, 4, 1), n))
# print(f(Permutation(4, 2, 5, 1, 3), n))
# print(f(Permutation(4, 5, 3, 1, 2), n))

# sigma = Permutation()
# print(sum(I_0[1]) - sum(apply_permutation(I_0, sigma)[1]))



f1 = None
for i in range(1, n+1):
    for j in I_0[1]:
        f1 = function_sum(f1, w_ij(i, j), b_coeff= j - i)

f1 = None
for i in range(1, n+1):
    for j in I_0[0]:
        f1 = function_sum(f1, w_ij(i, j), b_coeff= i - j)

# for part in [0, 1]:
#     for i in I_0[part]:
#         for j in I_0[part]:
#             print(i,j)
#             f1 = function_sum(f1, w_ij(i, j), b_coeff=j-i)

# for tau in Permutation.group(n):
#     if f(tau, n)-  f1(tau, n) != 0:
#         print(tau, f(tau, n), f1(tau, n))


# print(f1(Permutation() *  tau, n))
# print(f1(Permutation(1, 2, 4, 3, 5) * tau, n))
# print(f1(Permutation(1, 2, 5, 4, 3) * tau, n))
# print(f1(Permutation(1, 4, 3, 2, 5) * tau, n))
# print(f1(Permutation(1, 5, 3, 4, 2) * tau, n))
# print(f1(Permutation(1, 4, 5, 2, 3) * tau, n))
# print(f1(Permutation(4, 2, 3, 1, 5) * tau, n))
# print(f1(Permutation(5, 2, 3, 4, 1) * tau, n))
# print(f1(Permutation(4, 2, 5, 1, 3) * tau, n))
# print(f1(Permutation(4, 5, 3, 1, 2) * tau, n))

stat = distance_from_standard((3, 2), n)
# print(dft(f, n))

# print(dft(f1, n)[0])
# print(dft(f1, n)[1])
# print(dft(stat, n)[0])
# print(dft(stat, n)[1])

stab = [Permutation(1, 2, 3, 5, 4) *  sigma for sigma in Permutation.group(3)] + [sigma for sigma in Permutation.group(3)]
representatives = [Permutation(), Permutation(1, 2, 4, 3, 5), Permutation(1, 2, 5, 4, 3), Permutation(1, 4, 3, 2, 5), Permutation(1, 5, 3, 4, 2), Permutation(1, 4, 5, 2, 3), Permutation(4, 2, 3, 1, 5), Permutation(5, 2, 3, 4, 1), Permutation(4, 2, 5, 1, 3), Permutation(4, 5, 3, 1, 2)]

n = 5

# for rep in representatives:
#     avg = 0
#     for sigma in stab:
#         avg += stat(rep * sigma, n)/len(stab)
#     print(rep, avg, f1(rep, n))

# for sigma in Permutation.group(n):
#     print([sigma(i) for i in range(1, n+1)], " => ", stat(sigma, n))


# f = None
# n = 6
# for i in range(1, n//2+1):
#     for j in range(1, n//2+1):
#         if i != j:
#             for k in range(1, n//2+1):
#                 for l in range(1, n//2+1):
#                     if i > j and k < l:
#                         images = product([2*i, 2*i-1], [2*j, 2*j-1])
#                         indices = product([2*k, 2*k-1], [2*l, 2*l-1])
#                         ijkls = product(images, indices)
#                         for ij, kl in ijkls:
#                             f = function_sum(f, w_ij_kl(*ij, *kl))

# for m in dft(f, n):
#     print(m)
# print(f(Permutation(1, 5, 3, 4, 2, 6), n))
