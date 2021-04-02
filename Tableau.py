import numpy as np
from numpy.linalg import matrix_power
import math
import copy


class Tableau(object):
    def __init__(self, data):
        self.data = data
        self.size = 0
        self.shape = []
        for row in self.data:
            self.size += len(row)
            self.shape += [len(row)]

    def __repr__(self):
        s = ''
        for row in self.data:
            #s += (2*len(row) + 1) * '-' 
            #s += "\n"
            s += "|"
            for col in row:
                s+= str(col)
                s+="|"
            s += "\n"
        return s
    
    def compare_helper(self, other, n):
        for row_index in range(len(self.data)):
            if n in self.data[row_index] and n in other.data[row_index]:
                return self.compare_helper(other, n-1)
            elif n in self.data[row_index]:
                return True
            elif n in other.data[row_index]:
                return False

    def __gt__(self, other):
        return self.compare_helper(other, self.size)

    def __eq__(self, other):
        return self.data == other.data

    def find(self, k):
        for row in range(len(self.data)):
            for col in range(len(self.data[row])):
                if self.data[row][col] == k:
                    return row, col
        return "Not in Tableau"

    def signed_distance(self, k):
        k_row, k_col = self.find(k)
        content_k = k_col - k_row
        k_plus_one_row, k_plus_one_col = self.find(k+1) 
        content_k_plus_one = k_plus_one_col - k_plus_one_row
        return  content_k_plus_one - content_k

    def is_standard(self):
        for row in range(len(self.data)-1):
            for i in range(len(self.data[row+1])):
                if self.data[row][i] > self.data[row+1][i]:
                    return False
            for i in range(len(self.data[row])-1):
                if self.data[row][i] > self.data[row][i+1]:
                    return False
        for i in range(len(self.data[-1])-1):
            if self.data[-1][i] > self.data[-1][i+1]:
                    return False
        return True

    def switch(self, k):
        switched_tableau = Tableau(copy.deepcopy(self.data))
        k_row, k_col = self.find(k)
        k_plus_one_row, k_plus_one_col = self.find(k+1)
        switched_tableau.data[k_row][k_col] = k+1
        switched_tableau.data[k_plus_one_row][k_plus_one_col] = k
        return switched_tableau