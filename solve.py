from z3 import *



i = Int('i')
j = Int('j')
k = Int('k')

ind_one = [i, j, k]

q = Int('q')
r = Int('r')
s = Int('s')

ind_two = [q, r, s]

solv = Solver()

solv.add(i<j)
solv.add(j<k)
solv.add(q<r)
solv.add(r<s)

vals = [Or(x == 1, x == 2, x == 3, x == 4) for x in ind_one + ind_two]

for v in vals:
    solv.add(v)

solv.add(i == q)
solv.add(j == r)
solv.add(k != s)

print(solv.check())
print(solv.model())
