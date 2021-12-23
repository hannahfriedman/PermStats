## This file needs clean up ## 

from typing import *
from permutation import Permutation
from compress import tuples_w_relative_order
import math
import itertools

class Sij(object):
    def __init__(self, ind: tuple, img: tuple, parents=[]):
        self.ind = ind
        self.img = img

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

        sij = Sij(tuple(indices), tuple(images))
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
    def intersections(cls, sij: List['Sij'], num_to_intersect: int):
        if num_to_intersect == 1:
            return sij
        result = []
        sets_of_sets = itertools.combinations(sij, num_to_intersect)
        for s in sets_of_sets:
            new = s[0] + s[1]
            for i in range(2, len(s)):
                if type(new) != cls:
                    break
                new += s[i]
            if type(new) == Sij:
                result.append(new)
        return result

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
    def incl_excl(cls, sij: list, num_occur: int, add: bool, n: int):
        inter = cls.intersections(sij, num_occur)
        d = cls.summarize(inter)
        print(d)
        if d == {}:
            return 0
        if add:
            return cls.row_contr(d, n) + cls.incl_excl(sij, num_occur + 1, False, n)
        else:
            return -cls.row_contr(d, n) + cls.incl_excl(sij, num_occur + 1, True, n)

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
        x = 4
        d = {}
        c = cls.counts(sij)
        for s in c.keys():
            if len(s) == size:
                if c[s] == x:
                    print(s)
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

if __name__ == "__main__":
    sigma = Permutation(1,2,3)
    k = 3
    n = 6
    num_to_intersect = 2
    # sij = Sij.generate_all(sigma, k, n)
    sij5 = Sij.generate_all(sigma, k, 5)
    # sij6 = Sij.generate_all(sigma, k, 6)
    # sij7 = Sij.generate_all(sigma, k, 7)
    for i in range(2, 10):
        # print("n=6", i, Sij.summarize(Sij.intersections(sij6, i)))
        # print("n=7", i, Sij.summarize(Sij.intersections(sij7, i)))
        inter5 = Sij.intersections(sij5, i)
        # inter6 = Sij.intersections(sij6, i)
        summ5 =  Sij.summarize(inter5)
        # summ6 = Sij.summarize(inter6)
        print('n = 5', i, summ5)
        for size in summ5.keys():
            print('\t', size, ':', Sij.count_rel_order(inter5, size))
        # print('n = 6', i, summ6)
        # for size in summ6.keys():
        #     print('\t', size, ':', Sij.frequency(inter6, size))
