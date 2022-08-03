import numpy as np
from itertools import combinations
from permutation import Permutation
import matplotlib.pyplot as plt

from perm_stats import generate_w_ij
from random_walk import representation
from random_walk import w_ij_mat
from one_local_stats import make_subplots
from misc import function_to_vector
from misc import distribution

### Was originally created to determine sizes of orbits of one local functions under the action of Sn X Sn ###

class Equivalence(object):

    def __init__(self, funcs, mats):
        self.funcs = funcs
        self.reps = mats

    def append(self, new_funcs):
        self.funcs.append(new_funcs)

    def test_equivalent(self, mat):
        for m in self.reps:
            if (mat == m).all():
                return True
        return False

    def __len__(self):
        return len(self.funcs)


def wij_to_mat(rep, n):
    ''' returns a visual represenation of which entries in a matrix wij functions in list pick up'''
    result = np.zeros((n,n))
    for row in range(rep.shape[0]):
        for col in range(rep.shape[1]):
            result[row, col] = int(rep[row, col] >= 2)
    return result

n = 4
k = 6
equivalence_classes = []
for funcs in combinations(generate_w_ij(n), k):
    rep = sum([representation(wij, w_ij_mat, n, n, False) for wij in funcs])
    eq_cl_exists = False
    for eq_cl in equivalence_classes:
        if eq_cl.test_equivalent(rep):
            eq_cl_exists = True
            eq_cl.append(funcs)
            break
    if not eq_cl_exists:
        mats = []
        for sigma in Permutation.group(n):
            sigma_rep = w_ij_mat(sigma, n)
            for tau in Permutation.group(n):
                tau_rep = w_ij_mat(tau, n)
                new_mat = sigma_rep @ rep @ tau_rep
                new_mat_exists = False
                for m in mats:
                    if (m == new_mat).all():
                        new_mat_exists = True
                        break
                if not new_mat_exists:
                    mats.append(new_mat)
        equivalence_classes.append(Equivalence([funcs], mats))


for eq_cl in equivalence_classes:
    for i in range(len(eq_cl.reps)):
        print(eq_cl.reps[i])
        print(distribution(sum([np.array(function_to_vector(wij, n)) for wij in eq_cl.funcs[0]])))    
        plt.subplot(12, 12, i+1)
        plt.imshow(eq_cl.reps[i], vmin = 0, vmax = 16)
        plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False, labeltop = False)
    plt.show()


    
for eq_cl in equivalence_classes:
    print(len(eq_cl))
                
    
                          
