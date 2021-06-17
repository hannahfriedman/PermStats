import numpy as np
from permutation import Permutation
from typing import *
import copy
import misc
import math
from random_walk import rep
import perm_stats

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
        mat[indices.index(image), indices.index(index)] = 1
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

def w_rep(indices: List[int], images: List[int], n: int) -> np.array:
    dim = misc.falling_factorial(n, n-len(indices))
    mat = np.zeros((dim, dim))
    for sigma in Permutation.group(n):
        contained = True
        for i in range(len(indices)):
            if sigma(indices[i]) != images[i]:
                contained = False
                break
        if contained:
            mat += compress(sigma, len(indices), n)
    return mat

class Path:
    def __init__(self, *args):
        self.data = list(args)
    def __len__(self):
        return len(self.data)
    def __eq__(self, other):
        return same_contents(self.data, other.data)
    def __ge__(self, other):
        for i in other.data:
            if i not in self.data:
                return False
        return True
    def __le__(self, other):
        for i in self.data:
            if i not in other.data:
                return False
        return True
    def __gt__(self, other):
        return len(self) != len(other) and self >= other
    def __lt__(self, other):
        return len(self) != len(other) and self <= other
    def append(self, pair):
        self.data.append(pair)
    def __iter__(self):
        return iter(self.data)
    def __repr__(self):
        return str(self.data)

def num_excedances(mat: np.array) -> dict:
    d = {}
    all_sets = find_all_sets(mat.shape[0])
    for i in range(3):
        mat = np.rot90(mat)
    for row in range(mat.shape[0]):
        for col in range(mat.shape[1]):
            if row + col >= mat.shape[0]-1:
                mat[row, col] = 0
    print(mat)
    d[0] = 1
    # for i in range(mat.shape[0]-1):
        # visited = [(j, 0) for j in range(0, i)]
        # helper(mat, (mat.shape[0] - 2 - i,0), 0, d, [], [], mat.shape[0], all_sets, Path())
    helper(mat, d, find_all_sets(mat.shape[0]), 1, [])
    return d

def helper(mat: np.array, d: dict, all_sets: List[Path], k: int, path_tracker: List[Tuple[Path, int]]):
    d[k] = 0
    if k == mat.shape[0]-1:
        for path in all_sets:
            if len(path) == k:
                d[k] += 1
                path_tracker.append((path, d[k]))
    else:
        helper(mat, d, all_sets, k+1, path_tracker)
        print(f"k: {k}, path tracker: {path_tracker}")
        for path in all_sets:
            if len(path) == k:
                diff = 0
                for super_set in all_sets:
                    if path < super_set:
                        for pair in path_tracker:
                            if pair[0] == super_set:
                                diff += pair[1] 
                path_tracker.append((path, math.factorial(mat.shape[0] - k)-diff))
                if k == 1:
                    print(f"path: {path}, diff: {diff}, contribution: {math.factorial(mat.shape[0] - k) - diff}")
                d[k] += math.factorial(mat.shape[0] - k) - diff
                

def foo(mat: np.array, 
           index: Tuple[int, int], 
           depth: int,
           d: dict, 
           visited_rows: List[int],
           visited_cols: List[int],
           n: int, 
           all_sets: List[Path], 
           current_path: Path):
    amount_to_subtract = 0
    vr = copy.deepcopy(visited_rows)
    vc = copy.deepcopy(visited_cols)
    cp = copy.deepcopy(current_path)
    vr.append(index[0])
    vc.append(index[1])
    cp.append(index)
    print(f"Index: {index}, Depth: {depth+1}, Visited Rows: {visited_rows}, Visited Columns: {visited_cols}")
    for row in range(n-1):
        if row not in vr:
            for col in range(n-1):
                if col not in vc:
                    new = (row, col) 
                    if mat[new] != 0:
                        amount_to_subtract += helper(mat, new, depth+1, d, vr, vc, n, all_sets, cp)
            vr.append(row)
    factor = math.factorial(n - (depth +1))
    super_sets = [path for path in all_sets if path > cp]
    if len(super_sets) > 0:
        factor -= math.factorial(n - (depth +1) -1)
    if depth + 1 in d.keys():
        d[depth + 1] += factor
    else:
        d[depth + 1] = factor

    amount_to_subtract += (depth+1)*factor

    mat[index] -= amount_to_subtract
    return amount_to_subtract
    
def find_all_sets(n: int) -> List[Path]:
    if n == 2:
        return [Path((0,0))]
    existing_tuples = []
    for row in range(n-1):
        for col in range(n-1):
            if row + col == n-2:
                existing_tuples.append((row,col))
    result = find_all_sets(n-1) + [Path(pair) for pair in existing_tuples]
    for pair in existing_tuples:
        new_paths = []
        for path in result:
            add_pair = True
            for step in path:
                if step[0] == pair[0] or step[1] == pair[1]:
                    add_pair = False
                    break
            if add_pair:
                new_path = copy.deepcopy(path)
                new_path.append(pair)
                new_paths.append(new_path)
        result += new_paths
    return remove_duplicates(result)

def remove_duplicates(l: list):
    result = []
    for a in l:
        if a not in result:
            result.append(a)
    return result

def same_contents(a: list, b: list):
    if len(a) != len(b):
        return False
    for i in a:
        if i not in b:
            return False
    return True

if __name__ == "__main__":
    m = rep("exced", 7, False)
    print(num_excedances(m))
    d = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0}
    for sigma in Permutation.group(7):
        d[perm_stats.excedances(sigma, 7)] += 1
    print(d)
    # for bobo in find_all_sets(4):
    #     print(bobo)
 
    
""" 
([(1, 0), (0, 2)], 2) 
 ([(1, 1), (0, 2)], 2) 
 ([(2, 0), (0, 2)], 2)

1234
24 

1234
23

1234
2 4

1234
2134
"""