from permutation import Permutation
from sij import Sij
import numpy as np
from typing import *
from itertools import combinations

def all_ones(mat: np.array):
    for row in range(mat.shape[0]):
        for col in range(mat.shape[1]):
            if row == col and mat[row, col] != 0:
                return False
            elif row != col and mat[row, col] != 1:
                return False
    return True

def rec(mat: np.array, indices: list,  size: int) -> list:
    if size == 1:
        return indices
    prev = rec(mat, indicies, size - 1)
    



