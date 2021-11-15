from permutation import Permutation
import sys
from misc import choose
import matplotlib.pyplot as plt

n = int(sys.argv[1])
k = 3
count = 0
dist = []
for sigma in Permutation.group(n):
    if sigma.inversions() >= k:
        # print(sigma.inversions(), choose(sigma. inversions(), k))
        dist.append(sigma.inversions())
        count += choose(sigma.inversions(), k)
print(dist)
plt.boxplot(dist)
plt.show()
print(count)

