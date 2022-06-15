from dft import dft
from dft import dft_natural
from dft import matrix_rep_gen_natural
from dft import sort_partitions
from dft import generate_partitions
from small_dft import one_local_w_ij_kl_dft
from perm_stats import excedances
from perm_stats import length
from perm_stats import w_ij
from perm_stats import w_ij_kl
from tableau import Tableau

# i = 2
# j = 3
# k = 1
# l = 4
# n = 5
n = 5
for i in range(1, n+1):
    for j in range(1, n+1):
        if i != j:
            for k in range(1, n):
                for l in range(k + 1, n+1):
                    print(i, j, k, l)
                    if not (dft_natural(w_ij_kl(i, j, k, l), n)[3] == one_local_w_ij_kl_dft([(i, j, k, l)], n)[2]).all():
                        print(dft_natural(w_ij_kl(i, j, k, l), n)[3])
                        print(one_local_w_ij_kl_dft([(i, j, k, l)], n)[2])

n = 5
i = 3
partitions = sort_partitions(generate_partitions(n))
tableaux_by_shape = [Tableau.gen_by_shape(n, partition) for partition in partitions]
Tableau.sort(tableaux_by_shape[3])
print(tableaux_by_shape[i])

# print(dft_natural(length, n))
print(dft_natural(w_ij_kl(1, 3, 2, 1), n)[3])
# print(dft_natural(w_ij_kl(1, 4, 3, 1), n))
print(one_local_w_ij_kl_dft([(1, 3, 2, 1)], n)[2])

# for perm, mats in matrix_rep_gen_natural(n).items():
#    print(perm, mats[-3])
