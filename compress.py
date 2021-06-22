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
    Generates all tuples (lists) of length tuple_size from the set {1, ..., n} without repetition. 
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
    helper(mat.shape[0], d, find_all_sets(mat.shape[0]), 1, [])
    return d

def helper(n: int, d: dict, all_sets: List[Path], k: int, path_tracker: List[Tuple[Path, int]]):
    d[k] = 0
    if k == n-1:
        for path in all_sets:
            if len(path) == k:
                d[k] += 1
                path_tracker.append((path, d[k]))
    else:
        helper(n, d, all_sets, k+1, path_tracker)
        print(f"k: {k}, path tracker: {path_tracker}")
        for path in all_sets:
            if len(path) == k:
                diff = 0
                for super_set in all_sets:
                    if path < super_set:
                        for pair in path_tracker:
                            if pair[0] == super_set:
                                diff += pair[1] 
                path_tracker.append((path, math.factorial(n - k)-diff))
                if k == 1:
                    print(f"path: {path}, diff: {diff}, contribution: {math.factorial(n - k) - diff}")
                d[k] += math.factorial(n - k) - diff               

def num_patterns(sigma: Permutation, k: int, n: int) -> dict:
    col_indices = tuples_w_relative_order(range(1, k+1), k, n)
    row_indices = tuples_w_relative_order(sigma, k, n)
    max_path_length = 10
    all_sets = find_all_sets_patterns(sigma, col_indices, row_indices, k, n, max_path_length)
    d = {}
    path_tracker = []
    d[0] = 0
    patterns_helper(d, all_sets, 1, col_indices, path_tracker, max_path_length, n)
    d[0] = math.factorial(n) - sum(d.values())
    return d
    
def patterns_helper(d: dict, 
                    all_sets: List[Path], 
                    k: int, 
                    col_indices: List[List[int]],
                    path_tracker: List[Tuple[Path, int]], 
                    max_path_length: int, 
                    n: int):
    d[k] = 0
    if k == max_path_length:
        for path in all_sets:
            if len(path) == k:
                num_permutations = math.factorial(find_num_free(path, col_indices, n))
                d[k] += num_permutations
                path_tracker.append((path, num_permutations))
    else:
        patterns_helper(d, all_sets, k+1, col_indices, path_tracker, max_path_length, n)
        for path in all_sets:
            if len(path) == k:
                diff = 0
                num_permutations = math.factorial(find_num_free(path, col_indices, n))
                for super_set in all_sets:
                    if path < super_set:
                        for pair in path_tracker:
                            if pair[0] == super_set:
                                diff += pair[1]
                num_permutations -= diff
                d[k] += num_permutations
                path_tracker.append((path, num_permutations))

def find_num_free(path: Path, col_indices: List[List[int]], n: int) -> int:
    used_indices = []
    for pair in path:
        for i in col_indices[pair[1]]:
            if i not in used_indices:
                used_indices.append(i)
    return n - len(used_indices)

def find_all_sets_patterns(sigma: Permutation, 
                           col_indices: list, 
                           row_indices: list, 
                           k: int, n: int, 
                           max_path_length: int) -> List[Path]:
    pairs = []
    for i in col_indices:
        pairs += [Path((row_indices.index(j), col_indices.index(i))) for j in row_indices]
    if max_path_length == 1:
        return pairs 
    prev_sets = find_all_sets_patterns(sigma, col_indices, row_indices, k, n, max_path_length - 1)
    result = copy.deepcopy(prev_sets)
    for pair in pairs:
        for path in prev_sets:
            add_pair = True
            for step in path:
                # if step == (2, 9) and pair.data[0] == (6, 8):
                    # print(row_indices[step[0]])
                    # print(row_indices[pair.data[0][0]]) 
                    # print(col_indices[step[1]]) 
                    # print(col_indices[pair.data[0][1]])
                if not compatible(row_indices[step[0]], row_indices[pair.data[0][0]], col_indices[step[1]], col_indices[pair.data[0][1]]):
                    add_pair = False
                    break
            if add_pair:
                new_path = copy.deepcopy(path)
                new_path.append(pair.data[0])
                result.append(new_path)
        # result += new_paths
    return remove_duplicates(result)

def compatible(image1: list, image2: list, index1: list, index2: list) -> bool:
    if image1 == image2 or index1 == index2:
        return False
    for i in image1:
        if i in image2 and index1[image1.index(i)] != index2[image2.index(i)]:
            return False
    for i in index1:
        if i in index2 and image1[index1.index(i)] != image2[index2.index(i)]:
            return False
    return True

def tuples_w_relative_order(sigma: Union[Permutation, List[int]], k: int, n: int) -> List[List[int]]:
    possible_tuples = generate_tuples(k, n)
    permutation = sigma
    if type(sigma) == Permutation:
        permutation = permutation_to_list(sigma, k)
    return [i for i in possible_tuples if same_rel_order(permutation, i)]

def permutation_to_list(sigma: Permutation, n: int):
    return [sigma(i) for i in range(1, n+1)]

def same_rel_order(a: List[int], b: List[int]) -> bool:
    order_a = sorted(a)
    order_b = sorted(b)
    rel_order_a = [order_a.index(i) for i in a]
    rel_order_b = [order_b.index(i) for i in b]
    return rel_order_a == rel_order_b

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
    m = rep("exced", 5, False)
    # print(num_excedances(m))
    # print(m)
    # d = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0}
    # for sigma in Permutation.group(5):
    #     d[perm_stats.excedances(sigma, 5)] += 1
    # print(d)
    # for bobo in find_all_sets(4):
    #     print(bobo)
    # d = generate_tuples(3, 5)
    # for i in range(len(d)):
    #     d[i] = tuple(d[i])
    # print(dict(zip(d,range(len(d)))))
    # a = [1,2,3]
    # b = [2,3,1]
    k = 2
    n = 5
    sigma = Permutation(2, 1)
    print(num_patterns(sigma, k, n))
    # col_indices = tuples_w_relative_order(list(range(1, k+1)), k, n)
    # row_indices = tuples_w_relative_order(sigma, k, n)
    # print(col_indices)
    # print(row_indices)
    # print(compatible([1, 5, 4], [2, 4, 3], [1, 2, 5], [2, 4, 5]))
    # all_sets = find_all_sets_patterns(sigma, col_indices, row_indices, k, n, 1)
    # print(len(all_sets))
    # for path in all_sets:
    #     if len(path) == 2:
    #         index1 = col_indices[path.data[0][1]]
    #         image1 = row_indices[path.data[0][0]]
    #         index2 = col_indices[path.data[1][1]]
    #         image2 = row_indices[path.data[1][0]]
    #         # index3 = col_indices[path.data[2][1]]
    #         # image3 = row_indices[path.data[2][0]]
    #         print(compatible(image1, image2, index1, index2))
    # print(all_sets)
    # [(2, 9), (6, 8)]


 
    
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
     0          1          2          3          4          5          6          7          8          9
[[1, 3, 2], [1, 4, 2], [1, 4, 3], [1, 5, 2], [1, 5, 3], [1, 5, 4], [2, 4, 3], [2, 5, 3], [2, 5, 4], [3, 5, 4]]
[[1, 2, 3], [1, 2, 4], [1, 2, 5], [1, 3, 4], [1, 3, 5], [1, 4, 5], [2, 3, 4], [2, 3, 5], [2, 4, 5], [3, 4, 5]]

5!/2!3! = 10
10^2 = 100
  
"""