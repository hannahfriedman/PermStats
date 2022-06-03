import numpy as np
from dft import factor_sn
one_two = np.array([[-1, -1, -1], [0, 1, 0], [0, 0, 1]])
two_three = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
three_four = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])

# w34
print(one_two@two_three@three_four@one_two + two_three @ three_four + one_two @ two_three @ three_four + three_four + one_two @ three_four + two_three @ three_four @ one_two)

# w11
print(one_two@one_two + two_three + two_three @ three_four @ two_three + two_three @ three_four + three_four + three_four @ two_three)

# w12
print(one_two + one_two @ two_three + one_two@two_three@three_four@two_three + one_two @ two_three @ three_four + one_two @ three_four + one_two @ three_four @ two_three)

# w21
print(one_two + one_two @ three_four + two_three @ one_two + two_three @ three_four @ two_three @ one_two + three_four @ two_three @ one_two + two_three @ three_four @ one_two)

# (2, 2) partition
one_two = np.array([[1, 0], [-1, -1]])
two_three = np.array([[0, 1], [1, 0]])
three_four = np.array([[1, 0], [-1, -1]])

d = {(1,2): one_two, (2, 3): two_three, (3, 4): three_four}

# w34 23
print(two_three @ three_four + one_two @ two_three @ three_four)

def dft_wij_kl(perms):
    
    
