from itertools import combinations
from sage.combinat.permutation import Permutation

def is_increasing(sub_perm: list):
    for i in range(0, len(sub_perm)-1):
        if sub_perm[i] > sub_perm[i+1]:
            return False
    return True

def count_inc_seq(perm: Permutation, k: int, n: int):
    count = 0
    for tup in combinations(range(1, n+1), k):
        count += int(is_increasing([perm(i) for i in range(1, n+1) if i not in tup]))
    return count

def fac(n: int):
    if n == 0:
        return 1
    return n * fac(n - 1)

n = 5
d = {}

sigma = Permutation(list(range(1, n+1)))
for _ in range(fac(n)):
    # count = count_inc_seq(sigma, 2, n)
    count = sigma.longest_increasing_subsequence_length()
    if count in d.keys():
        d[count].append(sigma)
    else:
        d[count] = [sigma]
    sigma = sigma.next()

for key, val in d.items():
    print( key, [sigma.reduced_word_lexmin() for sigma in val])
    # print( key, len([sigma.to_cycles() for sigma in val]))
