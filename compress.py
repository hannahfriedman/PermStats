import numpy as np
from permutation import Permutation
from typing import *
import copy

def compress(sigma: Permutation, tuple_size: int, n: int) -> np.array:
    """
    Arguments: sigma - the permutation whose image should be represented as a matrix
               tuple_size - The length of the tuples indexing the matrix (the permutation representation 
                            is given by  tuple_size = 1, the regulst representatio is given by  tuple_size = n)
                n - n for which sigma is in Sn
    Returns a compression of sigma as a matrix. 
    """
    indices = generate_tuples(tuple_size, n)
    mat = np.zeros((len(indices), len(indices)))
    for index in indices:
        image = [sigma(i) for i in index]
        mat[indices.index(index), indices.index(image)] = 1
    return mat

def generate_tuples(tuple_size: int, n: int) -> List[List[int]]:
    """
    Generates all tuples of length tuple_size from the set {1, ..., n} without repetition. 
    """
    if tuple_size == 1:
        return [[i] for i in range(1, n+1)]
    else:
        shorter_tuples = generate_tuples(tuple_size-1, n)
        result = []
        for short_tuple in shorter_tuples:
            for i in range(1, n+1):
                if i not in short_tuple:
                    temp  = copy.deepcopy(short_tuple)
                    temp.append(i)
                    result.append(temp)
        return result


if __name__ == "__main__":
    sigma = Permutation(2,1)
    print(compress(sigma, 2, 3))

