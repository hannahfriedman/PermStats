from dft import projection
from dft import dft
from perm_stats import count_occurrences
from perm_stats import inc_seq_k
from perm_stats import average_tabloid_statistic
from perm_stats import fixed_points
from permutation import Permutation
from misc import function_to_vector
from typing import Callable
from typing import Tuple
import matplotlib.pyplot as plt

def projection_loss(f: Callable[[Permutation, int], int], n: int, *args) -> Tuple[float, float, float]:
    ''' Returns a tuple where the first number is the 2-norm squared of the function, 
    the second number is the 2-norm squared of the projection into spaces specified by  args,
    the third number is 2-norm squared of the remainder
    f is the function we are projecting
    f is defined in S_n
    args specifies which spaces to project into where 0 is the trivial, 1 is the (n - 1, 1), etc.
    '''
    function_size = sum([i**2 for i in function_to_vector(f, n)])
    projection_size = sum([i**2 for i in projection(f, n, *args)])
    return function_size, projection_size, function_size - projection_size

### Plot the ratio of the size of the projection to the size of the original function ###
n = 6
plt.style.use('ggplot')
plt.figure()
for k in range(n):
    # k determines which function we are projecting
    print(k)
    l = [] # Keep track of percentages for this function
    for i in range(11):
        # i determines who many spaces we are projecting into
        # Here the function is a tabloid statistic which is just a count of fixed k-tuples
        f, p, d = projection_loss(average_tabloid_statistic([n-k] + [1]*k, n),n, *list(range(i+1)))
        print(i, f, p, d, p/f)
        l.append(p/f)
    plt.gca().plot([i for i in range(11)], l)
plt.savefig('Figures/fix_k_mra.eps')
plt.show()

