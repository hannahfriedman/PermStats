from permutation import Permutation
import numpy as np
import cmath
from perm_stats import excedances
from perm_stats import length
from perm_stats import major_index
from perm_stats import w_ij

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

D8 = [Permutation(), Permutation(2,3,4,1), Permutation(3,4,1,2), Permutation(4,1,2,3), Permutation(2,1,4,3), Permutation(3,2,1,4), Permutation(4,3,2,1), Permutation(1,4,3,2)]


r = Permutation(2,3,4,5,6,1)
baby = [Permutation(), r, r*r, r*r*r, r*r*r*r, r*r*r*r*r]
D12 = baby +  [Permutation(2,1,6,5,4,3) * b for b in baby]

def distributions(f, n: int) -> list:
    if n == 4:
        H = [Permutation(), Permutation(2,3,4,1), Permutation(3,4,1,2), Permutation(4,1,2,3), Permutation(2,1,4,3), Permutation(3,2,1,4), Permutation(4,3,2,1), Permutation(1,4,3,2)]
    elif n == 5:
        baby = [Permutation(), Permutation(2,3,4,5,1), Permutation(3,4,5,1,2),  Permutation(4,5,1,2,3), Permutation(5,1,2,3,4)]
        H = baby  + [Permutation(1,5,4,3,2) * b for b in baby]
    elif n == 6:
        H = D12
    dude = {}
    for sigma in Permutation.group(n):
        shifted = [sigma.inverse() * h * sigma for h in H]
        stat = [f(h, n) for h in shifted]
        if str(stat) in dude.keys():
            dude[str(stat)].append(sigma)
        else:
            dude[str(stat)] = [sigma]
    return [eval(key) for key in dude.keys()]

    # print(dude)
    # for bruh in dude.keys():
    #     print(bruh)
    # print(len(dude))



#########################
## DFT for dihedral group ##

class Representation(object):
    def __init__(self, data: list):
        self.data = data
    def __mul__(self, other: "Representation"):
        l = []
        for i in range(len(self.data)):
            l.append(np.array(self.data[i] @ other.data[i]))
        return Representation(l)
    def __add__(self, other: "Representation"):
        l = []
        for i in range(len(self.data)):
            l.append(self.data[i] + other.data[i])
        return Representation(l)
    def __pow__(self, n: int):
#        print(self)
        if n == 1:
            return self
        else:
            return self * self**(n-1)
    def __str__(self):
        s = ''
        for m in self:
            s += str(m)
            s += '\n'
        return s
    def __iter__(self):
        return iter(self.data)
    def scale(self, n: int):
        new = [np.array(n*i) for i in self.data]
        return Representation(new)

# D8
r = Representation([np.array([[1]]), np.array([[1]]), np.array([[-1]]), np.array([[-1]]), np.array([[0, 1], [-1, 0]])])
s = Representation([np.array([[1]]), np.array([[-1]]), np.array([[1]]), np.array([[-1]]), np.array([[1, 0], [0, -1]])])
D8 = [r**4, r, r**2, r**3, s, s*r, s*r**2, s*r**3]

#D10
w = cmath.e**(2*cmath.pi*complex(0, 1)/5)
r = Representation([np.array([[1]]), np.array([[1]]), np.array([[w, 0],[0, w**(-1)]]), np.array([[w**2, 0],[0, w**(-2)]])])
s = Representation([np.array([[1]]), np.array([[-1]]), np.array([[0,1],[1,0]]), np.array([[0,1],[1,0]])])
D10 = [r**i for i in [5, 1, 2, 3, 4]] + [s*r**i for i in range(1, 6)]



n = 5
for dist in distributions(excedances, n):
    if n == 4:
        weighted = [D8[i].scale(dist[i]) for i in range(len(D8))]
    elif n == 5:
        weighted = [D10[i].scale(dist[i]) for i in range(len(D10))]
    DFT = weighted[0]
    for i in range(1, len(weighted)):
        DFT += weighted[i]
    DFT = [np.round_(a, 3) for a in DFT]
    print(Representation(DFT))


