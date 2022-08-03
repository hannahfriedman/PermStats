from perm_stats import w_I_J
from itertools import combinations
from misc import function_sum
from misc import function_to_vector
from dft import dft
from perm_stats import power_function
from perm_stats import fixed_points

n = 5
f = None

for k in range(2, 5):
    for I in combinations(list(range(1, n+1)), k):
        f = function_sum(f,w_I_J([I], [I]))
    print(function_to_vector(f, n))
    print(dft(f, n))

print(dft(power_function(fixed_points, k), n))
