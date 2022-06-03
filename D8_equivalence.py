import numpy as np
from itertools import combinations

n = 3


def index_to_matrix(index, n):
    m = np.zeros((n, n))
    for i in index:
        m[tuple(i)] = 1
    return m

def matrix_to_index(m):
    index = []
    for row in range(m.shape[0]):
        for col in range(m.shape[1]):
            if m[row, col] == 1:
                index.append((row, col))
    return index

def find_equivalent(mat):
    mats = [mat]
    mats.append(np.rot90(mats[-1]))
    mats.append(np.rot90(mats[-1]))
    mats.append(np.rot90(mats[-1]))
    mats.append(np.fliplr(mats[0]))
    mats.append(np.flipud(mats[0]))
    mats.append(np.transpose(mats[0]))
    mats.append(np.transpose(mats[2]))
    return [matrix_to_index(mat) for mat in mats]

class D8_Equivalence(object):
    def __init__(self, indices):
        self.indices = indices;
        self.mat = index_to_matrix(indices, n)
        self.equivalent = find_equivalent(self.mat)

    def __eq__(self, other):
        for blug in self.equivalent:
            if other.indices == blug:
                return True
        return False

    @classmethod
    def powerset_equivalence(cls, iterable):
        s = list(iterable)
        pow_set = []
        for i in range(0, len(s) + 1):
            for indices in combinations(s, i):
                # print(list(indices))
                if D8_Equivalence(list(indices)) not in pow_set:
                    pow_set.append(D8_Equivalence(list(indices)))
        print(len(pow_set))
        # for _ in pow_set:
        #     print(_.indices)
        print(pow_set[1].equivalent, pow_set[1].indices)
        print(pow_set[3].equivalent, pow_set[3].indices)        
        return pow_set




        
