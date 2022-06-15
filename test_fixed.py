from dft import dft
from dft import dft_natural
from perm_stats import fixed_points
from perm_stats import excedances
from perm_stats import length
from perm_stats import descent
from perm_stats import w_ij
from perm_stats import power_function

## Fixed points ## 
n = 5
d = {}
for k in range(61):
    f = power_function(fixed_points, k)
    triv = dft(f, n)[0][0][0]/120
    if triv in d.keys():
        d[triv].append(k)
    else:
        d[triv] = [k]
print(d)

n = 4
for i in range(1, n+1):
    for j in range(1, n+1):
        d = {}
        for k in range(13):
            f = power_function(w_ij(i, j), k)
            triv = dft(f, n)[0][0][0]
            if triv in d.keys():
                d[triv].append(k)
            else:
                d[triv] = [k]
        print(i,j,d)


# n = 5
# d = {}
# for k in range(61):
#     f = power_function(length, k)
# #     print(dft_natural(f, n))
#     triv = dft(f, n)[0][0][0]/20
#     if triv in d.keys():
#         d[triv].append(k)
#     else:
#         d[triv] = [k]
# print(d)

# n = 5
# d = {}
# for k in range(61):
#     f = power_function(descent, k)
# #     print(dft_natural(f, n))
#     triv = dft(f, n)[0][0][0]/4
#     if triv in d.keys():
#         d[triv].append(k)
#     else:
#         d[triv] = [k]
# print(d)

# n = 5
# d = {}
# for k in range(61):
#     f = power_function(excedances, k)
# #     print(dft_natural(f, n))
#     triv = dft(f, n)[0][0][0]/20
#     if triv in d.keys():
#         d[triv].append(k)
#     else:
#         d[triv] = [k]
# print(d)

    
