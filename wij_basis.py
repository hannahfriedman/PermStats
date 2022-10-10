from misc import function_sum_to_vector
from misc import function_to_vector
from misc import function_sum
from misc import matlab_syntax
from decompress import T_w_ij
from perm_stats import w_ij
from permutation import Permutation
import numpy as np



inc = [Permutation(1, 2, 3, 4), Permutation(1, 2, 4, 3), Permutation(1, 3, 2, 4),
       Permutation(1, 3, 4, 2), Permutation(1, 4, 2, 3),
       Permutation(2, 1, 3, 4), Permutation(2, 3, 1, 4),
       Permutation(2, 3, 4, 1), Permutation(3, 1, 2, 4), Permutation(4, 1, 2, 3)]

T = T_w_ij(4)

for sigma in inc:
    def f(tau, n):
        return sigma == tau
    print(np.reshape(T @ function_to_vector(f, 4), (4,4)))

for sigma in inc:
    def f(tau, n):
        return sigma == tau
    print(T @ function_to_vector(f, 4))
    

def create_matrix():
    m = np.zeros((24, 10))
    for ind in range(len(inc)):
        summands = [w_ij(inc[ind](i), i) for i in range(1, 5)]
        for i in range(1, 5):
            print(inc[ind](i), i)
        m[:, ind] = np.transpose(function_sum_to_vector(4, *summands)[0])
    return m

A = create_matrix()
B = np.linalg.inv(np.transpose(A) @ A) @ np.transpose(A) 

# print(B)

# print(np.linalg.inv(np.transpose(A) @ A))
# matlab_syntax(np.transpose(A) @ A)
# matlab_syntax(A)

# for i in range(1, 5):
#     for j in range(1, 5):
#         v = B @ function_to_vector(w_ij(i, j), 4)
#         print(i, j, v)
#         print(create_matrix() @ v)

# for sigma in Permutation.group(4):
    
