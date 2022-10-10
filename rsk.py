from tableau import Tableau
from permutation import Permutation
from perm_stats import w_ij_kl
from itertools import combinations

def rs(sigma: Permutation, n):
    tableau = []
    recording_tableau = []
    for i in range(1, n+1):
        # print(tableau, recording_tableau)
        tableau, recording_tableau = insert(tableau, recording_tableau, sigma(i), i)
    return Tableau(tableau), Tableau(recording_tableau)
        

def insert(tableau: list, recording_tableau: list, i: int, position: int):
    element_to_insert = i
    for row in range(len(tableau)):
        if tableau[row][-1] < element_to_insert:
            tableau[row].append(element_to_insert)
            recording_tableau[row].append(position)
            return tableau, recording_tableau
        else:
            for index in range(len(tableau[row])):
                if element_to_insert < tableau[row][index]:
                    # print(element_to_insert, tableau[row][index])
                    element_to_insert, tableau[row][index] = tableau[row][index], element_to_insert
                    break
    tableau.append([element_to_insert])
    recording_tableau.append([position])
    return tableau, recording_tableau

def is_increasing(sub_perm: list):
    for i in range(0, len(sub_perm)-1):
        if sub_perm[i] > sub_perm[i+1]:
            return False
    return True

def count_inc_seq(perm: Permutation, k: int, n: int):
    ''' return the  number of increasing subsequences of length n - k '''
    count = 0
    for tup in combinations(range(1, n+1), k):
#        print([perm(i) for i in range(1, n) if i not in tup])
        count += int(is_increasing([perm(i) for i in range(1, n+1) if i not in tup]))
    return count


def count_perms(n, i, j, k, l, m):
    perms = []
    for sigma in Permutation.group(n):
        if count_inc_seq(sigma, m, n) > 0:
            if w_ij_kl(i, j, k, l)(sigma, n):
                perms.append(sigma)
    return perms

n = 7
for sigma in Permutation.group(n):
    if count_inc_seq(sigma, 2, n) == n - 1:
        print(sigma)
# n = 5
# count = 0
# for i in range(1, n+1):
#     for j in range(1, n+1):
#         if i != j:
#             for k in range(1, n+1):
#                 for l in range(k, n+1):
#                     if k != l:
#                         # if len(count_perms(n, i, j, k, l, 2)) == 1:
#                         #     print(i, j, k, l)
#                         if i > j:
#                             if min(i - 2, k-1) + 1 + min(n-i, n- k- 1) < n - 2 and min(j-1, l - 2) + 1 + min(n - j - 1, n - l) < n - 2:
#                                 print(len(count_perms(n, i, j, k, l, 2)), i, j, k, l)
#                                 if len(count_perms(n, i, j, k, l, 2)) == 1:
#                                     count += 1
#                         else:
#                             if min(i-1, k-1) + 1 + min(n-i-1, n-k-1) < n - 2  and  min(j-2, l-2) + 1 + min(n-j, n-l) < n - 2 and min(i-1, k-1) + 1 + min(j - i -1, l - k -1) + 1 + min(n - j, n - l) < n - 2:# and len(count_perms(n, i, j, k, l, 2)) == 1:
#                                 print(len(count_perms(n, i, j, k, l, 2)), i, j, k, l)
#                                 if len(count_perms(n, i, j, k, l, 2)) == 1:
#                                     count += 1
                        #     if len(count_perms(n, i, j, k, l, 2)) > 1:
                        #         print(abs(i-k), abs(j-l))
                        #     if len(count_perms(n, i, j, k, l, 2)) == 1:
                        #         count += 1
                        #     print(i, j, k, l)
                            # print(*rs(count_perms(n, i, j, k, l, 2)[0], n))
                            # print(i, j, k, l)
#                            print(max(abs(i-k), abs(j-l)))
#                            print(abs(i+j - k-l))
  #                          print(abs(i - k) + abs(j - l))
#                            print(count_inc_seq(count_perms(n, i, j, k, l, 2)[0], 2, n))
#                        print(count_perms(n, i, j, k, l, 2))

#print(count)
i = 1
j = 2
m = 3
l = 4

n = 5
k = 2
# for sigma in Permutation.group(n):
    # P, Q = rs(sigma, n)
    # if P.shape[0] == n - k:
    #     print(count_inc_seq(sigma, k, n))
    #     print(P, Q)

#        print([sigma(i) for i in range(1, n+1)])
 
    # if w_ij_kl(i, j, m, l)(sigma, n):
    #     P, Q = rs(sigma, n)
    #     print([sigma(inx) for inx in range(1, n+1)])
    #     print(P, Q)
#         if P.shape[0] >= k:
#             print(count_inc_seq([sigma(i) for i in range(1, n+1)], n-k, n))
# #            print([sigma(i) for i in range(1, n+1)])
#             print(P, Q)

# print(rs(Permutation(4, 1, 2, 7, 6, 5, 8, 9, 3), 9))
