import numpy as np
import math
import copy
from typing import Type
from permutation import Permutation
from itertools import combinations

'''Tabloid class based on tableau class'''

class Tabloid(object):
    def __init__(self, data: list) -> None:
        self.data = [set(row) for row in data]
        self.size = 0
        self.shape = []
        for row in self.data:
            self.size += len(row)
            self.shape.append(len(row))

    def __repr__(self) -> str:
        s = '-' * self.shape[0] 
        for row in self.data:
            s += "\n"
            for val in row:
                s += str(val)
            s += '\n' + '-' * len(row)
        return s
                
                
    def __eq__(self, other):
        for row in range(len(self.data)):
            if self.data[row] != other.data[row]:
                return False
        return True

    def copy(self):
        return Tabloid([{i for i in row} for row in self.data])

    def apply_permutation(self, sigma: Permutation) -> 'Tabloid':
        return Tabloid([{sigma(i) for i in row} for row in self.data])

    def stabilizer(self, n: int):
        return [sigma for sigma in Permutation.group(n) if self.apply_permutation(sigma) == self]

    def descents(self):
        '''
        A descent occurs when i+1 is a in a row after i
        '''
        count = 0
        rows_checked = set()
        for row in self.data[:-1]:
            for i in row:
                # If i+1 is not the current row or in a previous row, it must be in one of the next rows
                if i+1 not in row and i+1 not in rows_checked:
                    count += 1
            rows_checked |= row
        return count

    def inversions(self):
        ''' An inversion occurs when j > i is in a row below i '''
        count = 0
        rows_checked = set()
        for row in self.data:
            for i in row:
                for j in range(i+1, self.size+1):
                    if j not in row and j not in rows_checked:
                        count += 1
            rows_checked |= row
        return count
        
    
    @staticmethod
    def gen_by_shape_helper(n_set: set, partition: list):
        if len(n_set) != sum(partition):
            raise Exception("Number of things do fill does not match shape")
        if len(partition) == 0:
            return []
        if len(partition) == 1:
            return [[n_set]]
        result = []
        for k_set in combinations(n_set, partition[0]):
            result += [[k_set] + t for t in Tabloid.gen_by_shape_helper(n_set - set(k_set), partition[1:])]
        return result


    @staticmethod
    def generate_by_shape(partition: list, n: int):
        return [Tabloid(t) for t in Tabloid.gen_by_shape_helper(set(range(1, n+1)), partition)]
