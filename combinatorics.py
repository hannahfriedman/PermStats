import numpy as np
from permutation import Permutation

def is_derangement(sigma, n):
    for i in range(1, n+1):
        if sigma(i) == i:
            return False
    return True

def generate_derangements(n):
    return [sigma for sigma in Permutation.group(n) if is_derangement(sigma, n)]

s4_derangements = generate_derangements(4)

n = 4
count = 0
for i in s4_derangements:
    for j in s4_derangements:
        if is_derangement(i * j, n):
            count += 1
print(count)

        
