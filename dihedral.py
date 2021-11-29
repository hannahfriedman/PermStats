from permutation import Permutation
from perm_stats import excedances
from perm_stats import length
from perm_stats import major_index
from perm_stats import w_ij

D8 = [Permutation(), Permutation(2,3,4,1), Permutation(3,4,1,2), Permutation(4,1,2,3), Permutation(2,1,4,3), Permutation(3,2,1,4), Permutation(4,3,2,1), Permutation(1,4,3,2)]

baby = [Permutation(), Permutation(2,3,4,5,1), Permutation(3,4,5,1,2),  Permutation(4,5,1,2,3), Permutation(5,1,2,3,4)]
D10 = baby  + [Permutation(1,5,4,3,2) * b for b in baby]

r = Permutation(2,3,4,5,6,1)
baby = [Permutation(), r, r*r, r*r*r, r*r*r*r, r*r*r*r*r]
D12 = baby +  [Permutation(2,1,6,5,4,3) * b for b in baby]

f = w_ij(2,2)
for n in [4,5,6]:
    if n == 4:
        H = D8
    elif n == 5:
        H = D10
    elif n == 6:
        H = D12
    dude = {}
    for sigma in Permutation.group(n):
        shifted = [sigma.inverse() * h * sigma for h in H]
        stat = [f(h, n) for h in shifted]
        if str(stat) in dude.keys():
            dude[str(stat)] += 1
        else:
            dude[str(stat)] = 1

    print(dude)
    # for bruh in dude.keys():
    #     print(bruh)
    print(len(dude))
