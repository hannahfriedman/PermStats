from decompress import TstarT_w_ij
from decompress import TstarT_w_ij_kl
from permutation import Permutation
from graphviz import Graph
from itertools import permutations
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
    return result

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
            

n = 5
mat = sum([rho_reg(sigma, n) for sigma in gen_rn(n)])
# mat = np.zeros((math.factorial(n), math.factorial(n)))
# for i in range(1, n):
#     print(gen_s_nk(n,i))
#     for sigma in gen_s_nk(n, i):
#         mat += rho_reg(sigma,n)
# for sigma in gen_s_nk(n, n-1):
#     mat += rho_reg(sigma, n)
print(mat)
print(TstarT_w_ij(n))
print(np.array_equal(mat, TstarT_w_ij(n)))

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
