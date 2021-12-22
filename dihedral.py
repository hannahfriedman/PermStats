###
## This file is for experimenting with permutation statistics for the dihedral group
###
from permutation import Permutation
import numpy as np
import cmath
from perm_stats import excedances
from perm_stats import length
from perm_stats import major_index
from perm_stats import w_ij

# Format numpy arrays when we print
np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

def distributions(f, n: int) -> list:
    '''
    Arguments: f, permutation statistic function in perm_stats.py
    n, the n of D2n
    Returns either a list of possible distributions according to the permutation statistic depending on the labeling of the vertices
    (Easily  modified to return a dictionary where the keys are distributions and the values are the number of labelings for which the distribution occurs)
    The distribution is the value of the permutation statistic on 1, r, r^2, ..., s, sr, sr^2, ... in that order.

    Note: this  function currently only works for n = 4, 5, 6, but should be easily modifiable
    '''
    # Create permutation representations
    if n == 4:
        H = [Permutation(), Permutation(2,3,4,1), Permutation(3,4,1,2), Permutation(4,1,2,3), Permutation(2,1,4,3), Permutation(3,2,1,4), Permutation(4,3,2,1), Permutation(1,4,3,2)]
    elif n == 5:
        baby = [Permutation(), Permutation(2,3,4,5,1), Permutation(3,4,5,1,2),  Permutation(4,5,1,2,3), Permutation(5,1,2,3,4)]
        H = baby  + [Permutation(1,5,4,3,2) * b for b in baby]
    elif n == 6:
        r = Permutation(2,3,4,5,6,1)
        baby = [Permutation(), r, r*r, r*r*r, r*r*r*r, r*r*r*r*r]
        H = baby +  [Permutation(2,1,6,5,4,3) * b for b in baby]

    # dude keeps track of how many times each distribution shows up for different labelings
    dude = {}
    # iterate over all possible labelings
    for sigma in Permutation.group(n):
        shifted = [sigma.inverse() * h * sigma for h in H] # All possible labellings are obtained by conjugating the original labeling
        stat = [f(h, n) for h in shifted] # Find the distribution for that labelling
        # Add distribution to the dictionary
        if str(stat) in dude.keys():
            dude[str(stat)].append(sigma)
        else:
            dude[str(stat)] = [sigma]
    return [eval(key) for key in dude.keys()]

############################
## DFT for dihedral group ##
############################

class Representation(object):
    '''
    This class is a list of matrices each of which corresponds to an irreducible representations
    It defines methods for add representations, exponentiating representations, and multiplying representations
    '''
    def __init__(self, data: list):
        '''
        The only attribute is the list of numpy arrays, each of which corresponds to an irreducible representation
        '''
        self.data = data
        
    def __mul__(self, other: "Representation"):
        '''
        Multiply representations pointwise
        '''
        l = []
        for i in range(len(self.data)):
            l.append(np.array(self.data[i] @ other.data[i]))
        return Representation(l)
    
    def __add__(self, other: "Representation"):
        '''
        Add representations pointwise
        '''
        l = []
        for i in range(len(self.data)):
            l.append(self.data[i] + other.data[i])
        return Representation(l)
    
    def __pow__(self, n: int):
        '''
        Recursively  exponentiate representations using mutliplication
        '''
        if n == 1:
            return self
        else:
            return self * self**(n-1)
        
    def __str__(self):
        '''
        Return string representation
        Each matrix is printed pretty with a line break between every matrix
        '''
        s = ''
        for m in self:
            s += str(m)
            s += '\n'
        return s
    
    def __iter__(self):
        return iter(self.data)

    def scale(self, n: int):
        '''
        Return a new representation that contains every matrix in self scaled by n
        '''
        new = [np.array(n*i) for i in self.data]
        return Representation(new)

def dihedral(n: int) -> list:
    '''
    Returns a representation  of the dihedral group of D2n
    '''
    # w is a primitive nth root of unity
    w = cmath.e**(2*cmath.pi*complex(0, 1)/n)
    # For odd n, these are the one dimensional representations of r and s and the number of two dimensional representations
    if n%2 == 1:
        onedr = [np.array([[1]]), np.array([[1]])]
        oneds = [np.array([[1]]), np.array([[-1]])]
        num_2d_reps = int((n-1)/2)
    # For even n, set the one dimensionsal representations of r and s and the number of two dimensional representations
    else:
        onedr += [np.array([[-1]]), np.array([[-1]])]
        oneds += [np.array([[1]]), np.array([[-1]])]
        num_2d_reps = int((n-2)/2)
    # Use w to generate the two dimensional representations for s and r and turn them into Representation objects
    r = Representation(onedr + [np.array([[w**k, 0],[0, w**(-k)]]) for k in range(1, num_2d_reps + 1)])
    s = Representation(oneds + num_2d_reps*[np.array([[0,1],[1,0]])])
    # Use r and s to generate the whole group
    D2n = [r**n] + [r**i for i in range(1, n)] + [s] + [s*r**i for i in range(1, n)]
    return D2n


if __name__ == "__main__":
    n = 5 # For now, n must be between 4 and 6
    f = w_ij(3,5)
    D2n = dihedral(n)
    # Print the DFT of f on D2n
    for dist in distributions(f, n):
        # Create a new list of represenations of D2n scaled by the distribution
        weighted = [D2n[i].scale(dist[i]) for i in range(len(D2n))]
        
        # Add up all the different scaled representations (can't use sum because empty sum is 0)
        DFT = weighted[0]
        for i in range(1, len(weighted)):
            DFT += weighted[i]

        DFT = [np.round_(a, 3) for a in DFT] # Makes the DFT easier to read
        print(dist)
        print(Representation(DFT))


