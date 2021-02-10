import numpy as np
import math


def matrix_rep(n):
    partitions = generate_partitions(n)
    tableaux_by_shape = [tableaux_shape(n, partition) for partition in partitions]
    rho = {}
    for i in range(1,n):
        representation = []
        for shape in tableaux_by_shape:
            sort_tableaux(n, shape)
            rep = np.zeros((len(shape), len(shape)))
            print(shape)
            for index in range(len(shape)):
                rep[index, index] = 1/(shape[index].signed_distance(i))
                switched = shape[index].switch(i)
#               if switched.is_standard():
 #                   switched_index = 0
  #                  for j in range(len(shape)):
   #                     if shape[j] == switched:
    #                        switched_index = j
     #                       break
      #              rep[switched_index][index] = math.sqrt(1 - (shape[index].signed_distance(i))^2)
            representation.append(rep)
        rho["(" + str(i) + "," + str(i+1) + ")"] = representation
    return rho


def generate_partitions(n):
    ans = []
    used = []
    if n == 1:
        return [[1]]
    elif n == 0:
        return [[]]
    for x in range(1, n):
        ans+=[[x] + part for part in generate_partitions(n-x)]
    return remove_dubs(ans) + [[n]]

def remove_dubs(partition):
    for part in partition:
        part.sort()
        part.reverse()
    return rem_dubs_helper(partition)

def rem_dubs_helper(partition):
    if len(partition) == 1:
        return partition
    for i in range(1, len(partition)):
        if partition[0] == partition[i]:
            return rem_dubs_helper(partition[1:])
    return [partition[0]] + rem_dubs_helper(partition[1:])



def generate_tableaux(n):
    if n == 1:
        return [Tableau([[1]])]
    else:
        prev_tab = generate_tableaux(n-1)
        ans = [Tableau(tab.data + [[n]]) for tab in prev_tab]
        for tab in prev_tab:
            for i in range(len(tab.data)):
                if tab.size == 1:
                    ans += [Tableau([tab.data[i] + [n]])]
                elif i == 0 or len(tab.data[i-1]) > len(tab.data[i]):
                    ans += [Tableau(tab.data[0:i] + [tab.data[i] + [n]] + tab.data[i+1:])]
        return ans



def tableaux_shape(n, partition):
    ans = []
    tableaux = generate_tableaux(n)
    for tab in tableaux:
        if len(tab.data) == len(partition):
            for row_index in range(len(tab.data)):
                matching = True
                if len(tab.data[row_index]) != partition[row_index]:
                    matching = False
            if matching:
                ans += [tab]
    return ans
                
def sort_tableaux(n, tableaux):
    switched = True
    while switched:
        switched = False
        for i in range(len(tableaux) - 1):
            if tableaux[i] < tableaux[i+1]:
                tableaux[i], tableaux[i+1] = tableaux[i+1], tableaux[i]
                switched = True
    return




class Tableau(object):
    def __init__(self, data):
        self.data = data
        self.size = 0
        for row in self.data:
            self.size += len(row)

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
        content_k = k_row - k_col
        k_plus_one_row, k_plus_one_col = self.find(k+1) 
        content_k_plus_one = k_plus_one_row - k_plus_one_col
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
        switched_tableau = Tableau(self.data)
        k_row, k_col = self.find(k)
        k_plus_one_row, k_plus_one_col = self.find(k+1)
        switched_tableau.data[k_row][k_col] = k+1
        switched_tableau.data[k_plus_one_row][k_plus_one_col] = k
        return switched_tableau
