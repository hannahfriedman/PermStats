import numpy as np
from dft import sort_partitions
from dft import generate_partitions
from tableau import Tableau
from math import factorial
one_two = np.array([[-1, -1, -1], [0, 1, 0], [0, 0, 1]])
two_three = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
three_four = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])

# # w34
# print(one_two@two_three@three_four@one_two + two_three @ three_four + one_two @ two_three @ three_four + three_four + one_two @ three_four + two_three @ three_four @ one_two)

# # w11
# print(one_two@one_two + two_three + two_three @ three_four @ two_three + two_three @ three_four + three_four + three_four @ two_three)

# # w12
# print(one_two + one_two @ two_three + one_two@two_three@three_four@two_three + one_two @ two_three @ three_four + one_two @ three_four + one_two @ three_four @ two_three)

# # w21
# print(one_two + one_two @ three_four + two_three @ one_two + two_three @ three_four @ two_three @ one_two + three_four @ two_three @ one_two + two_three @ three_four @ one_two)

# # (2, 2) partition
# one_two = np.array([[1, 0], [-1, -1]])
# two_three = np.array([[0, 1], [1, 0]])
# three_four = np.array([[1, 0], [-1, -1]])

# d = {(1,2): one_two, (2, 3): two_three, (3, 4): three_four}

# # w34 23
# print(two_three @ three_four + one_two @ two_three @ three_four)

def two_local_w_ij_kl_dft(wijkls, n):
    partitions = sort_partitions(generate_partitions(n))
    tableaux_by_shape = [Tableau.gen_by_shape(n, partition) for partition in partitions]
    for i in range(0, 3):
        Tableau.sort(tableaux_by_shape[i])
    n_min_one_one = sum([n_minus_one_one(*wijkl, n) for wijkl in wijkls])
    for wijkl in wijkls:
        if type(n_minus_two_one_one(*wijkl, tableaux_by_shape[3], n)) == type(None):
            print(wijkl)
    n_min_two_one_one = sum([n_minus_two_one_one(*wijkl, tableaux_by_shape[3], n) for wijkl in wijkls])
    return [np.array([[factorial(n - 2) * len(wijkls)]]), n_min_one_one, np.zeros((n*(n-3)//2, n*(n-3)//2)),  n_min_two_one_one]

def n_minus_one_one(i, j, k, l, n):
    result = np.zeros((n-1, n-1))
    if l == 1:
        k, l = l, k
        i, j = j, i
    if k == 1:
        if j != 1:
            result[j - 2, l - 2] = factorial(n - 2)
        if i != 1:
            result[i - 2, l - 2] = - factorial(n - 2)
        for col in range(result.shape[1]):
            if col != l - 2:
                for row in range(result.shape[0]):
                    if row == i - 2:
                        result[row, col] = - factorial(n-2)
                    elif row != j - 2:
                        result[row, col] = factorial(n - 3)
    else:
        for index in range(n-1):
            if index+2 != i and index+2 != j:
                result[index, k-2] = - factorial(n-3)
                result[index, l-2] = - factorial(n-3)
            if i != 1:
                result[i-2, k - 2] = factorial(n-2)
            if j != 1:
                result[j-2, l - 2] = factorial(n-2)
    return result


def n_minus_two_one_one(i, j, k, l, tableaux, n):
    result = np.zeros((len(tableaux), len(tableaux)))
    if i > j:
        i, j = j, i
        k, l = l, k
    if k == 1:
        if i != 1:
            for t in range(len(tableaux)):
                r,s = tableaux[t].data[1][0], tableaux[t].data[2][0]
                if r == l:
                    for row in range(len(tableaux)):
                        if tableaux[row].data[1][0] == i and tableaux[row].data[2][0] == j: # e_i_j
                            result[row, t] = factorial(n - 2)
                        elif tableaux[row].data[1][0] == i:
                            result[row, t] = -factorial(n - 3)
                        elif tableaux[row].data[1][0] == j:
                            result[row, t] = factorial(n - 3)
                        elif tableaux[row].data[2][0] == i:
                            result[row, t] = factorial(n - 3)
                        elif tableaux[row].data[2][0] == j:
                            result[row, t] = -factorial(n - 3) 
                elif s == l:
                    for row in range(len(tableaux)):
                        if tableaux[row].data[1][0] == i and tableaux[row].data[2][0] == j: # e_i_j
                            result[row, t] = - factorial(n - 2)
                        elif tableaux[row].data[1][0] == i:
                            result[row, t] = factorial(n - 3)
                        elif tableaux[row].data[1][0] == j:
                            result[row, t] = - factorial(n - 3)
                        elif tableaux[row].data[2][0] == i:
                            result[row, t] = - factorial(n - 3)
                        elif tableaux[row].data[2][0] == j:
                            result[row, t] = factorial(n - 3) 
            return result
        else: # i = 1
            for t in range(len(tableaux)):
                r,s = tableaux[t].data[1][0], tableaux[t].data[2][0]
                if r == l:
                    for row in range(len(tableaux)):
                        if tableaux[row].data[2][0] == j:
                            result[row, t] = - factorial(n - 3)
                        elif tableaux[row].data[1][0] == j:
                            result[row, t] = factorial(n - 3)
                elif s == l:
                    for row in range(len(tableaux)):
                        if tableaux[row].data[2][0] == j:
                            result[row, t] = factorial(n - 3)
                        elif tableaux[row].data[1][0] == j:
                            result[row, t] = - factorial(n - 3)
            return result
    elif l == 1:
        if i != 1:
            for t in range(len(tableaux)):
                r,s = tableaux[t].data[1][0], tableaux[t].data[2][0]
                if s == k:
                    for row in range(len(tableaux)):
                        if tableaux[row].data[1][0] == i and tableaux[row].data[2][0] == j: # e_i_j
                            result[row, t] = factorial(n - 2)
                        elif tableaux[row].data[1][0] == i:
                            result[row, t] = -factorial(n - 3)
                        elif tableaux[row].data[1][0] == j:
                            result[row, t] = factorial(n - 3)
                        elif tableaux[row].data[2][0] == i:
                            result[row, t] = factorial(n - 3)
                        elif tableaux[row].data[2][0] == j:
                            result[row, t] = -factorial(n - 3) 
                elif r == k:
                    for row in range(len(tableaux)):
                        if tableaux[row].data[1][0] == i and tableaux[row].data[2][0] == j: # e_i_j
                            result[row, t] = - factorial(n - 2)
                        elif tableaux[row].data[1][0] == i:
                            result[row, t] = factorial(n - 3)
                        elif tableaux[row].data[1][0] == j:
                            result[row, t] = - factorial(n - 3)
                        elif tableaux[row].data[2][0] == i:
                            result[row, t] = - factorial(n - 3)
                        elif tableaux[row].data[2][0] == j:
                            result[row, t] = factorial(n - 3) 
            return result
        else: # i = 1
            for t in range(len(tableaux)):
                r,s = tableaux[t].data[1][0], tableaux[t].data[2][0]
                if s == k:
                    for row in range(len(tableaux)):
                        if tableaux[row].data[2][0] == j:
                            result[row, t] = - factorial(n - 3)
                        elif tableaux[row].data[1][0] == j:
                            result[row, t] = factorial(n - 3)
                elif r == k:
                    for row in range(len(tableaux)):
                        if tableaux[row].data[2][0] == j:
                            result[row, t] = factorial(n - 3)
                        elif tableaux[row].data[1][0] == j:
                            result[row, t] = - factorial(n - 3)
            return result
    # k, l has a standard tableau associated with it
    for t in range(len(tableaux)):
        r,s = tableaux[t].data[1][0], tableaux[t].data[2][0]
        if r == k and s == l:
            if i != 1: # no need to check j, since i < j now
                for row in range(len(tableaux)):
                    if tableaux[row].data[1][0] == i and tableaux[row].data[2][0] == j: # e_i_j
                        result[row, t] = factorial(n - 2)
                    elif tableaux[row].data[1][0] == i: # e_i_q, so i < q
                        result[row, t] = -factorial(n - 3)
                    elif tableaux[row].data[2][0] == j: # e_q_j, so q < j
                        result[row, t] = -factorial(n - 3)
                    elif tableaux[row].data[2][0] == i: # e_q_i, so q < i < j
                        result[row, t] = factorial(n - 3)
                    elif tableaux[row].data[1][0] == j: # e_j_q, so i < j < q
                        result[row, t] = factorial(n - 3)
                return result
            else:
                for row in range(len(tableaux)):
                    if tableaux[row].data[1][0] == j:
                        result[row, t] = factorial(n - 3)
                    elif tableaux[row].data[2][0] == j:
                        result[row, t] = - factorial(n - 3)
                return result
        elif r == l and s == k:
            if i != 1: # no need to check j, since i < j now
                for row in range(len(tableaux)):
                    if tableaux[row].data[1][0] == i and tableaux[row].data[2][0] == j: # e_i_j
                        result[row, t] = -factorial(n - 2)
                    elif tableaux[row].data[1][0] == i: # e_i_q, so i < q
                        result[row, t] = factorial(n - 3)
                    elif tableaux[row].data[2][0] == j: # e_q_j, so q < j
                        result[row, t] = factorial(n - 3)
                    elif tableaux[row].data[2][0] == i: # e_q_i, so q < i < j
                        result[row, t] = -factorial(n - 3)
                    elif tableaux[row].data[1][0] == j: # e_j_q, so i < j < q
                        result[row, t] = -factorial(n - 3)
                return result
            else:
                for row in range(len(tableaux)):
                    if tableaux[row].data[1][0] == j:
                        result[row, t] = - factorial(n - 3)
                    elif tableaux[row].data[2][0] == j:
                        result[row, t] = factorial(n - 3)
                return result
                
                        
                            
                            

        
    
