from sage.combinat.permutation import Permutation
from sage.combinat.symmetric_group_representations import SymmetricGroupRepresentations
from math import factorial

def num_cycles(sigma, n):
    return len(sigma.cycle_type())

def num_k_cycles(k, leq = False):
    def f(sigma, n):
        count = 0
        for i in sigma.cycle_type():
            count += int(i == k)
            if leq:
                count += int(i < k)
        return count
    return f

def adj_cycles(pi, n):
    return sum(1 for cycle in pi.cycle_tuples()
               if cycle == tuple(range(cycle[0], cycle[0] + len(cycle))))

def max_adj_cycle(pi, n):
    adj_cycle_lengths = [len(cycle) for cycle in pi.cycle_tuples()
                if cycle == tuple(range(cycle[0], cycle[0] + len(cycle)))]
    if len(adj_cycle_lengths) == 0:
        return 0
    return max(adj_cycle_lengths)


def lis(pi, n):
    return pi.longest_increasing_subsequence_length()

def lis_smooth(pi, n):
    return pi.longest_increasing_subsequences_number() * pi.longest_increasing_subsequence_length()

def dft(f, n):
    sigma = Permutation(list(range(1, n+1)))
    mats = []
    orth = SymmetricGroupRepresentations(n, "orthogonal")
    for rep in orth:
        mats.append(f(sigma, n) * rep(sigma))
    for _ in range(factorial(n) - 1):
        sigma = sigma.next()
        val = f(sigma, n)
#        print(f(sigma, n), sigma)
        for i in range(len(orth)):
            mats[i] += val * orth[i](sigma)
    return mats

    # for partition in iterator:
    #     mats.append(sum())

# for m in dft(max_adj_cycle, 4):
#     print(m)

for m in dft(lis_smooth, 4):
    print(m)

# for m in dft(num_k_cycles(1), 4):
#     print(m)

# for m in dft(num_k_cycles(2), 4):
#     print(m)

# for m in dft(num_k_cycles(3), 4):
#     print(m)
    
# for m in dft(num_k_cycles(3, leq = True), 4):
#     print(m)
    
