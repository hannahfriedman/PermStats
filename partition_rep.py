from partition import Partition
from permutation import Permutation
from dft import sort_partitions
from dft import generate_partitions
from dft import dft
from perm_stats import fixed_points
import numpy as np

### Starting to investigate partition algebras ###

n = 5
partitions = sort_partitions(generate_partitions(n))

def partition_character(partition: tuple, sigma: Permutation, n: int) -> np.array:
    partition_basis = Partition.generate_by_shape(partition, n)
    return sum([1 for partition in partition_basis if partition == partition.apply_permutation(sigma)])

def conjugacy_class_reps(n: int) -> list:
    partitions = sort_partitions(generate_partitions(n))
    result = []
    for part in partitions:
        count = 1
        sigma = Permutation()
        for p in part:
            sigma *= Permutation.cycle(*list(range(count, count+p)))
            count += p
        result.append(sigma)
    return result

def char_S5(partition: list, sigma: Permutation):
    if partition == [5]:
        return 1
    if partition == [4, 1]:
        return fixed_points(sigma, 5) - 1
    if partition == [3, 2]:
        return round(np.trace(dft((lambda tau, n: int(tau == sigma)), 5)[2]))
    if partition == [3, 1, 1]:
        return round(np.trace(dft((lambda tau, n: int(tau == sigma)), 5)[3]))
    if partition == [2, 2, 1]:
        return round(np.trace(dft((lambda tau, n: int(tau == sigma)), 5)[4]))
    if partition == [2, 1, 1, 1]:
        return round(np.trace(dft((lambda tau, n: int(tau == sigma)), 5)[5]))
    return round(np.trace(dft((lambda tau, n: int(tau == sigma)), 5)[6]))

def conj_class_size_s5(sigma: Permutation):
    if sigma == Permutation():
        return 1
    if sigma == Permutation(2, 1):
        return 10
    if sigma == Permutation(2, 3, 1):
        return 20
    if sigma == Permutation(2, 1, 4, 3):
        return 15
    if sigma == Permutation(2, 3, 4, 1):
        return 30
    if sigma == Permutation(2, 3, 1, 5, 4):
        return 20
    if sigma == Permutation(2, 3, 4, 5, 1):
        return 24
    

for p1 in partitions:
#    if p1 in [[2, 1, 1, 1], [2, 2, 1], [3, 2]]:
    print('\n', p1)
    for p2 in partitions:
            # for sigma in conjugacy_class_reps(n):
            #     print(sigma, partition_character(p1, sigma, n), char_S5(p2, sigma), conj_class_size_s5(sigma))
        print(p2, sum([partition_character(p1, sigma, n)*char_S5(p2, sigma)*conj_class_size_s5(sigma) for sigma in conjugacy_class_reps(n)])/120)
