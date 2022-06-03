import numpy as np
import math
import copy
from numpy.linalg import matrix_power
from typing import Type
from permutation import Permutation

class Tableau(object):
    def __init__(self, data: list) -> None:
        """
        Tableau constructor
        takes in a list of lists or ints called data that contains the rows of the tableau
        """
        self.data = data
        self.size = 0     # Total number of entries in tableau
        self.shape = []   # Length of each row
        for row in self.data:
            self.size += len(row)
            self.shape += [len(row)]

    def __repr__(self) -> str:
        s = ''
        for row in self.data:
            s += "|"
            for col in row:
                s+= str(col)
                s+="|"
            s += "\n"
        return s
    
    def compare_helper(self, other: "Tableau", n: int) -> bool:
        '''
        Helper function for > operator
        '''
        for row_index in range(len(self.data)):
            if n in self.data[row_index] and n in other.data[row_index]: # If both tableau have n in the same row, find n-1
                return self.compare_helper(other, n-1)
            elif n in self.data[row_index]: # If it's in a higher row in self, self is greater
                return True
            elif n in other.data[row_index]: # If it's in a higher row in other, other is greater
                return False

    def __gt__(self, other):
        return self.compare_helper(other, self.size)

    def __eq__(self, other):
        return self.data == other.data

    def find(self, k: int) -> (int, int):
        for row in range(len(self.data)):
            for col in range(len(self.data[row])):
                if self.data[row][col] == k:
                    return row, col
        return -1, -1

    def signed_distance(self, k: int) -> int:
        """
        Return the signed distance from k to k+1
        """
        # Finding info for k
        k_row, k_col = self.find(k)
        content_k = k_col - k_row
        # Finding infor for k+1
        k_plus_one_row, k_plus_one_col = self.find(k+1) 
        content_k_plus_one = k_plus_one_col - k_plus_one_row
        # Return signed distance
        return  content_k_plus_one - content_k

    def is_standard(self) -> bool:
        """
        returns true if called on a standard tableau, false otherwise
        """
        for row in range(len(self.data)-1):
            # Columns should be strictly increasing
            for i in range(len(self.data[row+1])):
                if self.data[row][i] > self.data[row+1][i]:
                    return False
            # Rows should be strictly increasing
            for i in range(len(self.data[row])-1):
                if self.data[row][i] > self.data[row][i+1]:
                    return False
        # The last row should also be strictly decreasing
        for i in range(len(self.data[-1])-1):
            if self.data[-1][i] > self.data[-1][i+1]:
                    return False
        return True

    def switch(self, k: int) -> "Tableau":
        """
        Return a deep copy of the tableau with k and k+1 switched
        Requires k to be in the tableau
        """
        # Make a copy of the tableau
        switched_tableau = Tableau(copy.deepcopy(self.data))
        # Find k and k+1
        k_row, k_col = self.find(k)
        k_plus_one_row, k_plus_one_col = self.find(k+1)
        # Switch k and k+1
        switched_tableau.data[k_row][k_col] = k+1
        switched_tableau.data[k_plus_one_row][k_plus_one_col] = k
        return switched_tableau

    def dist_from_standard(self):
        ''' 
        Computes the distance of tableau from being standard, i.e. the inversions in a tableau
        '''
        dist = 0
        for row in range(len(self.data)):
            for col in range(len(self.data[row])):
                # Find inversions at self.data[row][col]
                for i in range(col + 1, len(self.data[row])):
                    if self.data[row][col] > self.data[row][i]:
                        dist += 1
                if row < len(self.data) - 1:# and len(self.data[row + 1]) > col:
                    for j in range(row+1, len(self.data)):
                        if len(self.data[j]) <= col:
                            break
                        elif  self.data[row][col] > self.data[j][col]:
                            dist += 1
        return dist

    def apply_permutation(self, sigma):
        data = []
        for row in self.data:
            data.append([sigma(i) for i in row])
        return Tableau(data)

    def permutation_difference(self, other):
        permutation  = self.size * [0]
        for row in range(len(self.data)):
            for col in range(len(self.data[row])):
                permutation[self.data[row][col] - 1] = other.data[row][col]
        return Permutation(*permutation)
                

    @staticmethod
    def perm_to_tableau(sigma, partition: tuple) -> "Tableau":
        data = []
        prev = 0
        for part in partition:
            row = []
            for entry in range(part):
                row.append(sigma(1 + prev + entry))
            data.append(row)
            prev += part
        return Tableau(data)

    @staticmethod
    def generate_all_lambda(lamb: tuple, n: int) -> list:
        '''
        Generates all tableuax (including non standard) of shape lambda
        '''
        return [perm_to_tableau(sigma) for sigma in Permutation.group(n)]

    @staticmethod
    def generate_all(n: int) -> list:
        '''
        Generates all tableaux of a given size n recursively and puts them in a list
        '''
        if n == 1:
            return [Tableau([[1]])]
        else:
            prev_tab = Tableau.generate_all(n-1)
            ans = []
            for tab in prev_tab:
                for i in range(len(tab.data)):
            # adds n to end of row if allowed
                    if i == 0 or len(tab.data[i]) < len(tab.data[i-1]):
                        ans += [Tableau(tab.data[0:i] + [tab.data[i] + [n]] + tab.data[i+1:])]
            # adds n as a new row
            ans += [Tableau(tab.data + [[n]]) for tab in prev_tab]
        return ans
    
    @staticmethod
    def gen_by_shape(n: int, partition: list) -> list:
        '''
        Generates all tableaux of size n, that fit a given partition
        '''
        ans = []
        # Generate all tableaux
        tableaux_list = Tableau.generate_all(n)
        for tab in tableaux_list:
            # Filter talbeaux by shape: only add tableaux whose shape match the correct partition
            if tab.shape == partition:
                ans += [tab]
        return ans
      
    @staticmethod
    def sort(tableaux: list) -> None:
        """
        sorts a list of tableau
        all of them should have the same shape
        """
        switched = True
        while switched:
            switched = False
            for i in range(len(tableaux) - 1):
                if tableaux[i] < tableaux[i+1]:
                    tableaux[i], tableaux[i+1] = tableaux[i+1], tableaux[i]
                    switched = True
        return

    
