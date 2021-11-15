import sys
from misc import choose
import math

n = int(sys.argv[1])

count = 0

def two_sets_gen(n: int) -> int:
    count = 0
    for i in range(1, n+1):
        for j in range(1, n+1):
            for k in range(i+1, n+1):
                for l in range(1, n+1):
                    if l != j:
                        count += 1
    return count

def three_sets_gen(n: int) -> int:
    count = 0
    for i in range(1, n+1):
        for j in range(1, n+1):
            for k in range(i+1, n+1):
                for l in range(1, n+1):
                    if l != j:
                        for r in range(k+1, n+1):
                            for s in range(1, n+1):
                                if s != l and s != j:
                                    count += 1
    return count

def four_sets_gen(n: int) -> int:
    count = 0
    for i in range(1, n+1):
        for j in range(1, n+1):
            for k in range(i+1, n+1):
                for l in range(1, n+1):
                    if l != j:
                        for r in range(k+1, n+1):
                            for s in range(1, n+1):
                                if s != l and s != j:
                                    for p in range(r+1, n+1):
                                        for q in range(1, n+1):
                                            if q != l and q != s and q != j:
                                                count += 1

    return count

def k_sets_gen(n: int, k: int) -> int:
    return math.factorial(k) * choose(n, k) ** 2


# print(two_sets_gen(n))
# print(three_sets_gen(n))
# print(four_sets_gen(n))


# Compute row 2 num sets
def two_sets(n: int) -> int:
    count = 0
    for j in range(1, n):
        for i in range(j+1, n+1):
            for k in range(j+1, n):
                for l in range(k+1, n+1):
                    if i != l:
                        count += 1
    return count

def three_sets(n: int) -> int:
    count = 0
    # Compute row 3 num sets
    for j in range(1, n):
        for i in range(j+1, n+1):
            for k in range(j+1, n):
                for l in range(k+1, n+1):
                    for r in range(k+1, n):
                        for s in range(r+1, n+1):
                            if i != l and s != l and i != s:
                                #                            print(j, '->', i, ',', k, '->', l, ',', r, '->', s)   
                                count += 1
    return count

def four_sets(n: int) -> int:
    count = 0
    # Compute row 4 num sets
    for j in range(1, n):
        for i in range(j+1, n+1):
            for k in range(j+1, n):
                for l in range(k+1, n+1):
                    for r in range(k+1, n):
                        for s in range(r+1, n+1):
                            for p in range(r+1, n):
                                for q in range(p+1, n+1):
                                    if i != l and s != l and i != s and q != i and q != l and q != s:
                                        count += 1
    return count

#print(two_sets(n))
#print('three sets', three_sets(n))
#print(four_sets(n))

def stirling(n: int, k: int) -> int:
    result = 0
    for i in range(k+1):
        result += (-1)**i * choose(k, i) * (k - i)**n
    return result/math.factorial(k)

####################################
# intersecting 2-set intersections #
####################################

def second_row(n: int) -> int:
    return choose(3, 1) *  stirling(n, n-3) + choose(4, 2) * stirling(n, n-4)/2

# def third_row(n: int) -> int:
#     return stirling(n, n-3) + 
    
#print(stirling(6, 3))
#print(two_sets(n))
#print(second_row(n))

##################
### Inversions ###
##################

def row_two(n: int) -> int:
    count = 0
    four = 0
    first = 0
    for i in range(1, n):
        for j in range(i+1, n+1):
            for l in range(1, n):
                for k in range(l+1, n+1):
                    for p in range(i, n):
                        if p == i:
                            r = k
                            for q in range(j+1, n+1):
                                for s in range(1, r):
                                    if s != l:
    #                                    print(i, j, '->', k, l, ', ', p, q, '->', r, s)
                                        first += 1
                                        count += 1
                        elif p == j:
                            r = l
                            for q in range(p+1, n+1):
                                for s in range(1, r):
                                    if s != l:
   #                                     print(i, j, '->', k, l, ', ', p, q, '->', r, s)
                                        count += 1
                        else:
                            for q in range(p+1, n+1):
                                if q == i:
                                    s = k
                                    for r in range(s+1, n+1):
                                        if r != k:
  #                                          print(i, j, '->', k, l, ', ', p, q, '->', r, s)
                                            count += 1
                                elif q == j:
                                    s = l
                                    for r in range(s+1, n+1):
                                        if r != k:
 #                                           print(i, j, '->', k, l, ', ', p, q, '->', r, s)
                                            count += 1
                                else:
                                    for s in range(1, n):
                                        for r in range(s+1, n+1):
                                            if s != l and r != k and s != k and r != l:
#                                                print(i, j, '->', k, l, ', ', p, q, '->', r, s)
                                                count += 1
                                                four += 1
                                        
                        # for q in range(j+1, n+1):
                        #     for s in range(1, n):
                        #         for r in range(s+1, n+1):
                        #             if ((k != r and i != p) or (k == r and i == p)) and ((l != s and j != q) or (l == s and j == q)):

    print(four)            #                 count += 1
    print(first)
    return count - four


print(row_two(n))

count = 0
for i in range(1, 4):
    for k in range(3, 6):
        count += (n - i - 1)*(n-i)*(k-1)*(k-2)
print(count)

count = 0
for i in range(2, 5):
    for k in range(2, 5):
        count += (n-i)*(i-1)*(k-1)*(n-k)
print(count)

count = 0
for i in range(1, n+1):
    for l in range(1, n+1):
        if choose(n-i, 2) * choose(n-i-2, 2) * (n-l) * (n-l-1) * (l-1) * (n-l-1) > 0:
            count += choose(n-i, 2) * choose(n-i-2, 2) * (n-l) * (n-l-1) * (l-1) * (n-l-1) 
print(count)
# print(k_sets_gen(n, 3))
# print(k_sets_gen(n, 4))
# print(k_sets_gen(n, 5))
# print(k_sets_gen(n, 6))
