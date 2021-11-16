from permutation import Permutation
from graphviz import Graph

def cayley(n: int, S: list) -> None:
    g = Graph()
    for sigma in Permutation.group(n):
        g.node(str(sigma), str(sigma))
    for s in S:
        for sigma in Permutation.group(n):
            g.edge(str(sigma), str(sigma*s))
    g.render('graphs/test', view=True)

def generate_sij(n: int, i: int, j: int) -> list:
    sij = []
    for sigma in Permutation.group(n-1):
        sij.append(Permutation.from_cycles((i, n)) *  sigma *   Permutation.from_cycles((n, j)))
    return sij

S = generate_sij(4, 1, 3)

print(S)
cayley(3, S)
