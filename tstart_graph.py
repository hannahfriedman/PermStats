from decompress import TstarT_w_ij
from decompress import TstarT_w_ij_kl
from permutation import Permutation
from graphviz import Graph
from itertools import permutations
from perm_stats import w_ij, w_ij_kl
from numpy.linalg import eig
from misc import wolfram_syntax
import numpy as np
import math

def gen_graph(n):
    mat = TstarT_w_ij(n)
    g = Graph()
    sn = [sigma for sigma in Permutation.group(n)]
#    print(sn)
    for sigma in sn:
        g.node(str(sigma), str(sigma))
    for i in range(len(sn)):
        for j in range(i, len(sn)):
            if mat[i, j] == 1:
                g.edge(str(sn[i]), str(sn[j]))
    g.render('graphs/test' + str(n), view=True)


def rho_reg(sigma, n) -> np.array:
    result = np.zeros((math.factorial(n), math.factorial(n)))
    sn = [perm for perm in Permutation.group(n)]
    for row in range(len(sn)):
        for col in range(len(sn)):
            if sn[row] == sigma*sn[col]:
                result[row, col] = 1
                break
    return np.transpose(result)

def gen_s_nk(n, k):
    perms = permutations(range(1, n+1), k)
    real = []
    for p in perms:
        sort = sorted(p)
        mapping = {sort[i]: p[i] for i in range(len(p))}
        real.append([mapping[i] if i in p else i for i in range(1, n+1)])
    return [Permutation(*p) for p in real]

def gen_rn(n):
    result = []
    for i in range(1, n):
        result += [Permutation.cycle(i, n) * sigma * Permutation.cycle(i, n) for sigma in Permutation.group(n-1)]
    result += [sigma for sigma in Permutation.group(n-1)]
    return result
            

n = 3
#mat = sum([rho_reg(sigma, n) for sigma in gen_rn(n)])
### WIJ VERSION ###
mat = sum([rho_reg(sigma, n) for sigma in Permutation.group(n) if w_ij(1,1)(sigma, n) == 1])
for i in range(2, n+1):
    j = i
    mat += sum([rho_reg(sigma, n) for sigma in Permutation.group(n) if w_ij(i, j)(sigma, n) == 1])

### WIJKL VERSION ###
# mat = sum([rho_reg(sigma, n) for sigma in Permutation.group(n) if w_ij_kl(1,2,1,2)(sigma, n) == 1])
# for i in range(1, n+1):
#     for j in range(i+1, n+1):
#         if j != 2:
#             mat += sum([rho_reg(sigma, n) for sigma in Permutation.group(n) if w_ij_kl(i, j, i, j)(sigma, n) == 1])

print(mat)
print(TstarT_w_ij(n))
print(np.array_equal(mat, TstarT_w_ij(n)))
# print(np.round_(eig(TstarT_w_ij(n))[0], decimals=1))
# print(np.round_(eig(TstarT_w_ij(n))[1], decimals=1))
# print(TstarT_w_ij(n) @ np.array([[1], [1], [1], [1], [1], [1]]))

### OLD ###
# mat = np.zeros((math.factorial(n), math.factorial(n)))
# for i in range(1, n):
#     print(gen_s_nk(n,i))
#     for sigma in gen_s_nk(n, i):
#         mat += rho_reg(sigma,n)
# for sigma in gen_s_nk(n, n-1):
#     mat += rho_reg(sigma, n)

#gen_graph(4)


# n = 5
# sn = [sigma for sigma in Permutation.group(n)]
# mat = TstarT_w_ij(n)
# entries = []
# #for sigma in [Permutation(2,3,1), Permutation(3,1,2)]:
# #for sigma in [Permutation(2,1), Permutation(1,3,2), Permutation(3,2,1)]:
# for i in range(len(sn)):
#     for j in range(i, len(sn)):
#         zero = True
#         for sigma in Permutation.group(n-2):
#             if sn[i] == sigma *  sn[j]:
#                 zero = False
#                 break
#         if zero:
#             entries.append(mat[i,j])
#         #if sn[i] == sigma *  sn[j]:
#             #entries.append(mat[i, j])
# print(entries)
