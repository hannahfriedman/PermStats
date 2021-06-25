from typing import *
from compress import *
from graphviz import Digraph
from permutation import Permutation

def graph(all_sets: List[Path]) -> None:
    g = Digraph(comment="Inversions in S4")
    for path in range(len(all_sets)):
        g.node(str(all_sets[path]), str(all_sets[path]))
        for sub_set in range(path):
            if all_sets[path] > all_sets[sub_set]:
                g.edge(str(all_sets[path]), str(all_sets[sub_set]))
    g.render('test-output/poset', view=True) 

if __name__ == "__main__":
    k = 2
    n = 4
    max_path_length = 6
    sigma = Permutation(2, 1)
    col_indices = tuples_w_relative_order(range(1, k+1), k, n)
    row_indices = tuples_w_relative_order(sigma, k, n)
    all_sets = find_all_sets_patterns(sigma, col_indices, row_indices, k, n, max_path_length)
    graph(all_sets)