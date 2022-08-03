from dft import dft
from dft import inverse_dft
from dft import dft_matrix
from dft import dft_natural
from perm_stats import fixed_points
from perm_stats import excedances
from perm_stats import length
from perm_stats import descent
from perm_stats import w_ij
from perm_stats import w_I_J
from perm_stats import power_function
from perm_stats import count_cycles
from itertools import combinations
from misc import function_sum
from misc import function_to_vector
from permutation import Permutation

print('hi')
n = 6
k = 3
one_n = list(range(1, n+1))
f = None
for i, j, l in combinations(one_n, k):
    # f = function_sum(f, w_I_J([(j,), (l,), (i,)], [(i,), (j,), (l,)]),b_coeff=k)
    # f = function_sum(f, w_I_J([(l,), (i,), (j,)], [(i,), (j,), (l,)]),b_coeff=k)
    f = function_sum(f, w_I_J([{i, j, l}], [{i, j, l}]), b_coeff = 1)
    f = function_sum(f, w_I_J([{i, j}, {l}], [{i, j}, {l}]), b_coeff = -1)
    f = function_sum(f, w_I_J([{i, l}, {j}], [{i, l}, {j}]), b_coeff = -1)
    f = function_sum(f, w_I_J([{l, j}, {i}], [{l, j}, {i}]), b_coeff = -1)
for i in one_n:
    for j in one_n:
        for l in one_n:
            if i != j and j != l and i != l:
                f = function_sum(f, w_I_J([{l}, {j}, {i}], [{l}, {j}, {i}]), b_coeff = 1/3)
# for i in one_n:
#     f = function_sum(f, w_ij(i, i))

for sigma in Permutation.group(n):
    # if f(sigma, n) - power_function(fixed_points, k)(sigma, n) != 0:
    #     print(sigma, f(sigma, n), power_function(fixed_points, k)(sigma, n))
    if round(f(sigma, n) - count_cycles(k)(sigma, n),5) != 0:
        print(sigma, f(sigma, n), count_cycles(k)(sigma, n))    

space = {6: ['(n)', '(n - 1, 1)', '(n - 2, 2)', '(n - 2, 1, 1)', '(n - 3, 3)', '(n - 3, 2, 1)', '(n - 3, 1, 1, 1)', '(n - 4, 2, 2)', '(n-4, 2, 1, 1)', '(n-4, 1, 1,1 , 1)', '(1, 1, 1, 1, 1, 1)'], 7: ['(n)', '(n - 1, 1)', '(n - 2, 2)', '(n - 2, 1, 1)', '(n - 3, 3)', '(n - 3, 2, 1)', '(n - 3, 1, 1, 1)', '(n-4, 3, 1)', '(n - 4, 2, 2)', '(n-4, 2, 1, 1)', '(n-4, 1, 1,1 , 1)', '(n-5, 2, 2, 1)', '(n-5, 2, 1, 1, 1)', '(n-5, 1, 1, 1, 1, 1)', '(1, 1, 1, 1, 1, 1, 1)']}
k = 3
n = 5
a = 1
fix = dft(fixed_points, n)
fix_k = dft(power_function(fixed_points, k), n)
d = {}
for j in range(2, k+1):
    if k%j == 0:
        d[j] = dft(count_cycles(j), n)


for i in range(len(fix)):
    print('----------------------------------------')
    if n in d.keys():
        if i < len(space[n]):
            print(space[n][i])
    current = fix[i]**a
    print('dim', current.shape[0])
    print(1,current[0,0])
    for j in range(2, k+1):
        if j in d.keys():
            print(j, (j**a * d[j][i]**a)[0,0])
            if j == 10:
                current += -j**a * d[j][i]**a
            else:
                current += j**a * d[j][i]**a
    print(current[0,0])
    print((fix_k[i]**a)[0,0])
print('----------------------------------------')


n = 6
k = 2
fix_k = dft(power_function(fixed_points, k), n)
# print(fix_k)
# k = 3
# for n in range(3, 9):
#     count_cycle = dft(count_cycles(k), n)
#     print(k, [m[0,0] for m in count_cycle])
# fix = dft(fixed_points, n)
# print(fix_k)
# print(count_cycle)
# print(fix)

# for i in range(len(fix)):
#     print(fix_k[i]**2)
#     print(k**2 * count_cycle[i]**2 + fix[i]**2)    

# print([str(sigma) for sigma in Permutation.group(n)])
# v = inverse_dft([m**2 for m in dft(power_function(fixed_points, 2), n)], n)
# print(function_to_vector(fixed_points, n))

# counter= 0
# for sigma in Permutation.group(n):
#     print([sigma(i) for i in range(1,n+1)], '=>', int(v[counter]))
#     counter += 1
