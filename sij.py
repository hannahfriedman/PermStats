from typing import *
import statistics
import numpy as np
from permutation import Permutation
from compress import tuples_w_relative_order
from misc import choose
import math
import sys
import matplotlib.pyplot as plt
import networkx as nx

patt_freq = []

class Sij(object):
    def __init__(self, ind: tuple, img: tuple, parents=[]):
        self.ind = ind
        self.img = img
        self.parents = parents

    def __getitem__(self, index: int) -> int:
        try:
            return self.img[self.ind.index(index)]
        except Exception:
            raise Exception(f"{ind}'s mapping is not specified in {self}.")

    def __repr__(self):
        return str(self.ind) + " -> " + str(self.img)

    def __eq__(self, other) -> bool:
        return self.ind == other.ind and self.img == other.img

    def __len__(self):
        return len(self.ind)

    def __hash__(self):
        return hash(repr(self))
    
    def same_parents(self, other) -> bool:
        if self != other or len(self.parents) != len(other.parents) or self.parents == [] or other.parents == []:
            return False
        for i in self.parents:
            if i not in other.parents:
                return False
        return True

    def unique_parents(self, sij: list) -> bool:
        for i in sij:
            if self == i and self.same_parents(i):
                return False
        return True
    
    def __add__(self, other) -> "Sij":
        """
        Returns intersection of Sij sets
        """
        indices = []
        images = []
        parents = []
        for i in self.ind:
            if (i in other.ind and self[i] != other[i]) or (self[i] in other.img and i != other.ind[other.img.index(self[i])]):
                return None
            indices.append(i)
            images.append(self[i])
        for j in other.ind:
            if j not in indices:
                indices.append(j)
                images.append(other[j])

        if self.parents == []:
            parents.append(self)
        else:
            parents += self.parents
        if other.parents == []:
            parents.append(other)
        else:
            parents += other.parents
        sij = Sij(tuple(indices), tuple(images), parents=parents)
        # print(sij)
        sij.reorder()
        return sij

    def reorder(self):
        indices = [i for i in self.ind]
        indices.sort()
        images = [self[i] for i in indices]
        self.ind = tuple(indices)
        self.img = tuple(images)

    @classmethod
    def generate_all(cls, sigma: Permutation, k: int, n: int) -> List["Sij"]:
        indices = tuples_w_relative_order(Permutation(*range(1, k+1)), k, n)
        images = tuples_w_relative_order(sigma, k, n)
        result = []
        for i in indices:
            for j in images:
                result.append(Sij(tuple(i), tuple(j)))
        return result

    @classmethod
    def generate_fancy_graph(cls, sij: List['Sij']) -> nx.graph:
        graph = nx.Graph()
        for i in sij:
            graph.add_node(i)
            for j in graph:
                if j != i:
                    new = i + j
                    if type(new) == Sij:
                        graph.add_edge(i, j)
        return graph
        
    @classmethod
    def generate_graph(cls, sij: List['Sij']):
        graph = {}
        for i in sij:
            graph[i] = []
            for j in graph.keys():
                if j != i:
                    new = i + j
                    if type(new) == Sij:
                        graph[j].append(i)
                        graph[i].append(j)
        return graph

    @classmethod
    def generate_graph_mat(cls, sij: List['Sij']):
        graph = np.zeros((len(sij), len(sij)))
        for i in range(len(sij)):
            for j in range(i, len(sij)):
                if i != j:
                    new = sij[i] + sij[j]
                    if type(new) == Sij:
                        graph[i, j] = 1
                        graph[j, i] = 1
        return graph


    @classmethod
    def intersections(cls, sij: List['Sij'], num_to_intersect: int, memory=True, previous=None):
        if num_to_intersect == 1:
            return sij
        result = []
        num_able_to_intersect = {i: 0 for i in sij}
        if num_to_intersect == 2:
            for i in range(len(sij)):
                for j in range(i+1, len(sij)):
                    new = sij[i] + sij[j]
                    if type(new) == Sij:
                        if len(new) ==  5:
                            num_able_to_intersect[sij[i]] += 1
                            num_able_to_intersect[sij[j]] += 1
#                        print(sij[i], "+", sij[j], "=", sij[i] + sij[j])
                        add = False
                        if not memory:#  and new not in result:
                            if new not in result:
                                new.parents = []
                                result.append(new)
                        else:
                            result.append(new)
                            # new.parents = []
                            # add = True
                        # if memory or add:
                        #     result.append(new)
            v = num_able_to_intersect.values()
            print(v)
            print(statistics.mean(v))
            print(statistics.median(v))
            plt.boxplot(v, meanline=True)
            plt.show()
        else:
            if type(previous) == list:
                prev = previous
            else:
                prev = cls.intersections(sij, num_to_intersect-1, memory=True)
            d = cls.summarize(prev)
            for size in d.keys():
                print(size, cls.count_rel_order(prev, size))
                patt_freq.append(cls.count_rel_order(prev, size))
            for i in sij:
                for j in prev:
                    new = i+j
                    if type(new) == Sij and i not in j.parents:
                        add = False
                        if not memory:
                            if new not in result:
                                new.parents = []
                                result.append(new)
                        elif new.unique_parents(result):
                            result.append(new)
                        # if not memory and new not in result:
                        #     new.parents = []
                        #     add = True
                        # if (not memory and add) or new.unique_parents(result):
                        #     result.append(new)
        return result

    @classmethod
    def print_intersections(cls, sij: List['Sij'], num_to_intersect: int) -> None:
        data = Sij.intersections(sij, num_to_intersect)
        # print(data)
        for i in data:
            print(i, ":", i.parents)

    @classmethod
    def summarize(cls, sij):
        result = {}
        for i in sij:
            if len(i) in result.keys():
                result[len(i)] += 1
            else:
                result[len(i)] = 1
        return result

    @classmethod
    def incl_excl(cls, sij: list, num_occur: int, add: bool, n: int, previous=None):
        inter = cls.intersections(sij, num_occur, previous=previous)
        d = cls.summarize(inter)
        # print(d, inter)
        print(d)
        # print(Sij.counts(inter))
        s = Sij.counts(inter)
        m = {}
        for i in s:
            if len(i) == 6:
                if s[i] in m.keys():
                    m[s[i]] += 1
                else:
                    m[s[i]] = 1
        print(m)
            # print(i, ":", s[i])
        # print([s for s in inter  if len(s) == 4])
        # cls.print_intersections(sij, num_occur)
        if d == {}:
            return 0
        if add:
            return cls.row_contr(d, n) + cls.incl_excl(sij, num_occur + 1, False, n, previous=inter)
        else:
            return -cls.row_contr(d, n) + cls.incl_excl(sij, num_occur + 1, True, n, previous=inter)

    @classmethod
    def row_contr(cls, d: dict, n: int) -> int:
        result = 0
        # print(d)
        for l in d.keys():
            result += d[l] * math.factorial(n - l)
        # print(result)
        return result
        
    @classmethod
    def at_least(cls, sigma, num_occ, k, n):
        sij = cls.generate_all(sigma, k, n)
        if num_occ != 1:
            sij = cls.intersections(sij, num_occ, memory=False)
        return cls.incl_excl(sij, 1, True, n)    

    @classmethod
    def counts(cls, sij: list) -> dict:
        d = {}
        for s in sij:
            if s in d.keys():
                d[s] += 1
            else:
                d[s] = 1
        return d

    @classmethod
    def frequency(cls, sij: list, size: int) -> dict:
        d = {}
        c = cls.counts(sij)
        for s in c.keys():
            if len(s) == size:
                if c[s] in d.keys():
                    d[c[s]] += 1
                else:
                    d[c[s]] = 1
        return d


    @classmethod
    def count_rel_order(cls, sij: list, size: int) -> dict:
        result = {}
        for pattern in sij:
            if len(pattern) == size:
                rel_order = reduce(list(pattern.img))
                if rel_order in result.keys():
                    result[rel_order] += 1
                else:
                    result[rel_order] = 1
        return result


def reduce(pattern: list) -> tuple:
    order = sorted(pattern)
    return tuple([order.index(i) for i in pattern])

def consolodate(patt_freq: List[dict]) -> dict:
    d = {}
    for i in patt_freq:
        for key,  value in i.items():
            if key in d.keys():
                d[key].append(value)
            else:
                d[key] = [value]
    return d

def sum_diff(values: List[int]) -> int:
    values = [(-1)**i * values[i] for i in range(len(values))]
    return sum(values)

def count_complete_graphs(graph: dict, size: int, previous = None):
    if size == 1:
        lonely_sij = []
        for key, value in graph.items():
            if value == []:
                lonely_sij.append([key])
        return lonely_sij
    if previous == None:
        previous = count_complete_graphs(graph, size-1)
    k_size = []
    for sub_graph in previous:
        valid = True
        for node in sub_graph:
            for candidate in graph[node]:
                for other_node in sub_graph:
                    if candidate not in graph[other_node]:
                        valid = False
                        break
                if valid:
                    k_size.append(sub_graph + [candidate])
    return k_size

def count_complete_all(graph: dict, n: int) -> dict:
    result = {1:0}
    node_counts = {key:math.factorial(n - len(key)) for key in graph.keys()}
    for node, count in node_counts.items():
        if count > 0:
            groups = [[] for x in range(count)]
            for neighbor0 in graph[node]:
#                node_counts[neighbor0] -= 1
                for group in groups:
                    valid = True
                    for neighbor in group:
                        if neighbor not in graph[neighbor0]:
                            valid = False
                            break
                    if valid:
                        group.append(neighbor0)
                        break
            total_not_one = 0
            for group in groups:
                for neighbor in group:
                    graph[neighbor] = [x for x in graph[neighbor] if (x != node and x not in group)]
                if group != []:
                    s = node
                    for neighbor in group:
                        s += neighbor
                    if len(group)+1 in result.keys():
                        result[len(group)+1] += math.factorial(n - len(s))
                    else:
                        result[len(group)+1] = math.factorial(n - len(s))
                    node_counts[neighbor] -= math.factorial(n - len(s))
                    total_not_one += math.factorial(n - len(s))
            result[1] += count - total_not_one
    return result

if __name__ == "__main__":
    if len(sys.argv) > 1:
        args = [int(i) for i in sys.argv[1:]]
        sigma = Permutation(*args)
        k = len(args)
    else:
        sigma = Permutation(2,3,1)
        k = 3
    n = 5
    sij = Sij.generate_all(sigma, k, n)
    graph = Sij.generate_fancy_graph(sij)
#    print(graph.edges())
    # for c in nx.algorithms.clique.find_cliques(graph):
    #     print(len(c))
#    print(graph)





















#    print(count_complete_all(graph, n))
# print(Sij.at_least(sigma, 1, k, n))
    # d = {}
    # for key, value in consolodate(patt_freq).items():
    #     print(key, ':', value, '=', sum_diff(value)/choose(n, len(key))**2)
    #     if len(key) in d.keys():
    #         d[len(key)].append(sum_diff(value)/choose(n, len(key))**2)
    #     else:
    #         d[len(key)] = [sum_diff(value)/choose(n, len(key))**2]
    # for key, value in d.items():
    #     d[key] = sorted(value)
    # print(d)
    # inter = Sij.intersections(sij, num_to_intersect, memory=False)
    # print(Sij.summarize(inter))
    # Sij.print_intersections(inter, 2)
    # for i in range(1, 10):
        # print(Sij.at_least(sigma, i, k, n))

    # for i in range(6, 7):
        # print(Sij.at_least(sigma, 1, k, i))
    # for i in range(5, 7):
        # sij = Sij.generate_all(sigma, k, i)
        # inter = Sij.intersections(sij, num_to_intersect)
        # inter = Sij.intersections(sij, 2)
        # for j in range(3, 20): 
            # print('row', j-1, 'n = ', i, ":", Sij.summarize(inter))
            # inter = Sij.intersections(sij, j, previous=inter)
        # inter = Sij.intersections(sij, num_to_intersect)
        # summ = Sij.summarize(inter)
        # print(i, ":", summ)
        # for size in summ.keys():
        #     print('\t', size, ':', Sij.frequency(inter, size))
        
# print(Sij.incl_excl(sij, i, True, n))

