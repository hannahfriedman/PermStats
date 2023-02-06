from compress import *
from rsk import count_inc_seq
from random_walk import rep_w_ij_kl
from perm_stats import w_ij_kl
import numpy as np
from permutation import Permutation
from itertools import permutations
from math import factorial

def perm_w_fewest(val:list, pos:list, n: int) -> Permutation:
    comp = [i for i in range(1, n+1) if i not in val]
    perms = [tup for tup in permutations(comp, len(comp))]
    perms.reverse()
    for perm in perms:
        sigma = Permutation(*merge(val, perm, pos, n))
        if count_inc_seq(sigma, 2, n) > 0:
            return sigma
        

def merge(val: list, comp: list, pos: list, n: int):
    result = []
    curr_pos = 0
    curr_comp = 0
    for i in range(1, n+1):
        if  curr_pos < len(pos):
            print(i + 1 == pos[curr_pos])
        if curr_pos < len(pos) and len(result) + 1 == pos[curr_pos]:
            result.append(val[curr_pos])
            curr_pos += 1
        else:
            result.append(comp[curr_comp])
            curr_comp += 1
    return result

def fill_in(val:list, pos:list, n: int) -> Permutation:
    order = []
    curr_pos = 0
    for i in range(1, n+1):
        while curr_pos < len(pos) and len(order) + 1 == pos[curr_pos]:
            order.append(val[curr_pos])
            curr_pos += 1
        if i not in val:
            order.append(i)
    return Permutation(*order)

def find_coeff(c, d, e, f, n):
    M = rep_w_ij_kl(c, d, e, f, n)
    zero = np.zeros(M.shape)
    indices = [tup for tup in permutations(range(1, n+1), 2)]
    print(indices)
    while not (M == zero).all():
        mag = np.absolute(M)
        a, b = np.unravel_index(np.argmax(mag), mag.shape)
        if M[a,b] == mag[a,b]:
            coeff = -1
        else:
            coeff = 1
        print(a, b)
        i, j = indices[a]
        k, l = indices[b]
        if l < k:
            k, l = l, k
            i, j = j, i
        print([i, j], [k, l])
        sigma = perm_w_fewest([i, j], [k, l], n)
        M = M  + coeff * compress(sigma, 2, n)
        print(i, j, k, l)
        print(M)
    print('Yay')
    



def inc_seq_in_supp(i, j, k, l, n):
    perms = []
    for sigma in Permutation.group(n):
        if sigma(k) == i and sigma(l) == j and count_inc_seq(sigma, 2, n) > 0:
            perms.append(sigma)
    return perms


d = {Permutation(): -1.0, Permutation(2, 1): 1.0, Permutation(1, 3, 2): 2.0, Permutation(3, 1, 2): -0.0, Permutation(2, 3, 1): -1.0, Permutation(3, 2, 1): 1.0, Permutation(1, 2, 4, 3): 1.0, Permutation(2, 1, 4, 3): -1.0, Permutation(1, 4, 2, 3): -0.0, Permutation(4, 1, 2, 3): 2.0, Permutation(2, 4, 1, 3): 1.0, Permutation(4, 2, 1, 3): -1.0, Permutation(1, 3, 4, 2): -1.0, Permutation(3, 1, 4, 2): 1.0, Permutation(1, 4, 3, 2): 1.0, Permutation(4, 1, 3, 2): -1.0, Permutation(3, 4, 1, 2): -1.0, Permutation(4, 3, 1, 2): 1.0, Permutation(2, 3, 4, 1): 1.0, Permutation(3, 2, 4, 1): -1.0, Permutation(2, 4, 3, 1): -1.0, Permutation(4, 2, 3, 1): 1.0, Permutation(3, 4, 2, 1): 2.0, Permutation(1, 2, 3, 5, 4): 0.0, Permutation(2, 1, 3, 5, 4): -0.0, Permutation(1, 3, 2, 5, 4): -0.0, Permutation(3, 1, 2, 5, 4): 0.0, Permutation(2, 3, 1, 5, 4): 0.0, Permutation(1, 2, 5, 3, 4): -0.0, Permutation(2, 1, 5, 3, 4): 0.0, Permutation(1, 5, 2, 3, 4): 0.0, Permutation(5, 1, 2, 3, 4): 0.0, Permutation(2, 5, 1, 3, 4): -0.0, Permutation(5, 2, 1, 3, 4): -0.0, Permutation(1, 3, 5, 2, 4): 0.0, Permutation(3, 1, 5, 2, 4): -0.0, Permutation(1, 5, 3, 2, 4): -0.0, Permutation(5, 1, 3, 2, 4): -0.0, Permutation(3, 5, 1, 2, 4): -0.0, Permutation(5, 3, 1, 2, 4): 0.0, Permutation(2, 3, 5, 1, 4): -0.0, Permutation(2, 5, 3, 1, 4): 0.0, Permutation(5, 2, 3, 1, 4): 0.0, Permutation(1, 2, 4, 5, 3): -0.0, Permutation(2, 1, 4, 5, 3): 0.0, Permutation(1, 4, 2, 5, 3): 0.0, Permutation(4, 1, 2, 5, 3): -0.0, Permutation(2, 4, 1, 5, 3): -0.0, Permutation(1, 2, 5, 4, 3): 0.0, Permutation(1, 5, 2, 4, 3): -0.0, Permutation(5, 1, 2, 4, 3): -0.0, Permutation(1, 4, 5, 2, 3): -0.0, Permutation(4, 1, 5, 2, 3): 0.0, Permutation(1, 5, 4, 2, 3): 0.0, Permutation(5, 1, 4, 2, 3): 0.0, Permutation(4, 5, 1, 2, 3): 0.0, Permutation(5, 4, 1, 2, 3): 0.0, Permutation(2, 4, 5, 1, 3): -0.0, Permutation(1, 3, 4, 5, 2): 0.0, Permutation(3, 1, 4, 5, 2): -0.0, Permutation(1, 4, 3, 5, 2): -0.0, Permutation(4, 1, 3, 5, 2): 0.0, Permutation(3, 4, 1, 5, 2): -0.0, Permutation(1, 3, 5, 4, 2): -0.0, Permutation(1, 5, 3, 4, 2): 0.0, Permutation(5, 1, 3, 4, 2): 0.0, Permutation(1, 4, 5, 3, 2): 0.0, Permutation(3, 4, 5, 1, 2): 0.0, Permutation(2, 3, 4, 5, 1): 0.0, Permutation(3, 2, 4, 5, 1): -0.0, Permutation(2, 4, 3, 5, 1): -0.0, Permutation(4, 2, 3, 5, 1): 0.0, Permutation(3, 4, 2, 5, 1): 0.0, Permutation(2, 3, 5, 4, 1): -0.0, Permutation(2, 5, 3, 4, 1): 0.0, Permutation(5, 2, 3, 4, 1): 0.0, Permutation(2, 4, 5, 3, 1): 0.0, Permutation(3, 4, 5, 2, 1): 0.0}

def find_coeff_2(c, d, e, f, n, max_iter = 3):
    M = rep_w_ij_kl(c, d, e, f, n)
    zero = np.zeros(M.shape)
    indices = [tup for tup in permutations(range(1, n+1), 2)]
    which_perm = np.zeros(M.shape)
    supports = {}
    for i in range(1, n+1):
        for j in range(1, n+1):
            if i != j:
                for k in range(1, n+1):
                    for l in range(1, n+1):
                        supports[(i,j,k,l)] = inc_seq_in_supp(i, j, k, l, n)
                        
    for key, perms in supports.items():
        if len(perms) > 1:
            for sigma in perms:
                for small_indices, small_perms in supports.items():
                    if len(small_perms) ==1 and sigma in small_perms:
                        #   print(perms, small_perms)
                        perms.remove(sigma)
                        break
            # if (c,d,e,f) in [(1, 2, 4, 5), (1, 4, 4, 5), (2, 1, 4, 5), (2, 5, 3, 5), (3, 1, 2, 4), (3, 1, 2, 5), (3, 2, 2, 4), (3, 4, 1, 4), (3, 4, 1, 5), (3, 5, 2, 3), (4, 1, 3, 4), (4, 1, 3, 5), (4, 1, 4, 5), (4, 2, 3, 4), (4, 2, 3, 5), (4, 3, 1, 2), (4, 5, 1, 3), (5, 1, 2, 4), (5, 1, 2, 5), (5, 1, 3, 5), (5, 2, 3, 5)]:
            #     perms = [perms[-1]] + perms[:-1]
            # perms.reverse()
    
    num_inc_seq_all = np.zeros(M.shape)
    for alpha in range(num_inc_seq_all.shape[0]):
        for beta in range(num_inc_seq_all.shape[1]):
            num_inc_seq_all[alpha, beta] = len(supports[(*indices[alpha], *indices[beta])])
#    print(num_inc_seq_all)

    linear_combo = {sigma:0 for sigma in Permutation.group(n)}
    
    counter = 0
    # print(M)
    while not (M == zero).all():
        for num_inc_seq in range(1, 8):
            for count in range(1, num_inc_seq + 1):
            # for count in range(1, num_inc_seq+1):
                count = num_inc_seq - count + 1
                # print(count)
                for a in range(M.shape[0]):
                    for b in range(M.shape[1]):
                        if M[a,b] != 0 and num_inc_seq_all[a,b] == count:
                            # print(a, b)
                            i, j = indices[a]
                            k, l = indices[b]
                            perms = supports[(i,j,k,l)]
                            sigma = perms[int(which_perm[a,b])]
                            which_perm[a,b] = (int(which_perm[a,b]) + 1)%len(perms)
                            linear_combo[sigma] += M[a,b]
                            M = M - M[a,b] * compress(sigma, 2, n)
                            if linear_combo == d:
                                print("HOLY FUCKING SHIT")
                # for key, val in linear_combo.items():
                #     print(val, key)
                # print()

 #       print('------------------------------------------------------------------------------------------------------------------------------------------------------------------')
        counter += 1
        # print(i, j, k, l)

        if counter == max_iter:
            # for m in M:
            #      print(list(m))
            print(c, d, e, f)
            print(':(')
            return 1
    return 0
  #  print('Yay')


def find_coeff_3(c, d, e, f, n, max_iter = 3):
    M = rep_w_ij_kl(c, d, e, f, n)
    zero = np.zeros(M.shape)
    indices = [tup for tup in permutations(range(1, n+1), 2)]
    which_perm = np.zeros(M.shape)
    supports = {}
    for i in range(1, n+1):
        for j in range(1, n+1):
            if i != j:
                for k in range(1, n+1):
                    for l in range(1, n+1):
                        supports[(i,j,k,l)] = inc_seq_in_supp(i, j, k, l, n)
                        
    for key, perms in supports.items():
        if len(perms) > 1:
            for sigma in perms:
                for small_indices, small_perms in supports.items():
                    if len(small_perms) ==1 and sigma in small_perms:
                        #   print(perms, small_perms)
                        perms.remove(sigma)
                        break
#            perms.reverse()
    
    num_inc_seq_all = np.zeros(M.shape)
    for alpha in range(num_inc_seq_all.shape[0]):
        for beta in range(num_inc_seq_all.shape[1]):
            num_inc_seq_all[alpha, beta] = len(supports[(*indices[alpha], *indices[beta])])
    print(num_inc_seq_all)
#    print(num_inc_seq_all)

    linear_combo = {}
    
    counter = 0
    # print(M)
    while not (M == zero).all():

        for num_inc_seq in range(1, 8):
            for count in range(1, num_inc_seq + 1):
            # for count in range(1, num_inc_seq+1):
                count = num_inc_seq - count + 1
                # print(count)
                for a in range(M.shape[0]):
                    for b in range(M.shape[1]):
                        if M[a,b] != 0 and num_inc_seq_all[a,b] == count:
                            # print(a, b)
                            i, j = indices[a]
                            k, l = indices[b]
                            perms = supports[(i,j,k,l)]
                            sigma = perms[int(which_perm[a,b])]
                            which_perm[a,b] = (int(which_perm[a,b]) + 1)%len(perms)
                            if sigma in linear_combo.keys():
                                linear_combo[sigma] += M[a,b]
                            else:
                                linear_combo[sigma] = M[a,b]
                            M = M - M[a,b] * compress(sigma, 2, n)

        mag = np.absolute(M)
        a, b = np.unravel_index(np.argmax(mag), mag.shape)
        if M[a,b] == mag[a,b]:
            coeff = -1
        else:
            coeff = 1
        i, j = indices[a]
        k, l = indices[b]
        perms = supports[(i,j,k,l)]
        sigma = perms[int(which_perm[a,b])]
        which_perm[a,b] = (int(which_perm[a,b]) + 1)%len(perms)
        # if sigma in linear_combo.keys():
        #     linear_combo[sigma] += M[a,b]
        # else:
        #     linear_combo[sigma] = M[a,b]
        M = M - np.sign(M[a,b]) * compress(sigma, 2, n)

        # count = 1
        # for a in range(M.shape[0]):
        #     for b in range(M.shape[1]):
        #         if M[a,b] != 0 and num_inc_seq_all[a,b] == count:
        #             i, j = indices[a]
        #             k, l = indices[b]
        #             perms = supports[(i,j,k,l)]
        #             sigma = perms[int(which_perm[a,b])]
        #             which_perm[a,b] = (int(which_perm[a,b]) + 1)%len(perms)
        #             if sigma in linear_combo.keys():
        #                 linear_combo[sigma] += M[a,b]
        #             else:
        #                 linear_combo[sigma] = M[a,b]
        #             M = M - M[a,b] * compress(sigma, 2, n)

                    
        counter += 1
        if counter == max_iter:
            for m in M:
                print(list(m))
            print(c, d, e, f)
            print(':(')
            return 1
    return 0




def find_coeff_4(c, d, e, f, n, max_iter = 1):
    M = rep_w_ij_kl(c, d, e, f, n)
    zero = np.zeros(M.shape)
    indices = [tup for tup in permutations(range(1, n+1), 2)]
    which_perm = np.zeros(M.shape)
    supports = {}
    for i in range(1, n+1):
        for j in range(1, n+1):
            if i != j:
                for k in range(1, n+1):
                    for l in range(1, n+1):
                        supports[(i,j,k,l)] = inc_seq_in_supp(i, j, k, l, n)
                        
    for key, perms in supports.items():
        if len(perms) > 1:
            for sigma in perms:
                for small_indices, small_perms in supports.items():
                    if len(small_perms) ==1 and sigma in small_perms:
                        perms.remove(sigma)
                        break
    
    num_inc_seq_all = np.zeros(M.shape)
    for alpha in range(num_inc_seq_all.shape[0]):
        for beta in range(num_inc_seq_all.shape[1]):
            num_inc_seq_all[alpha, beta] = len(supports[(*indices[alpha], *indices[beta])])
#    print(num_inc_seq_all)

    linear_combo = {sigma:0 for sigma in Permutation.group(n)}
    
    counter = 0
    # print(M)
    while not (M == zero).all():
        #for num_inc_seq in range(1, 8):
        #for count in range(1, num_inc_seq + 1):
        # for count in range(1, num_inc_seq+1):
        # count = num_inc_seq - count + 1
                # print(count)
        for a in range(M.shape[0]):
            for b in range(M.shape[1]):
                if M[a,b] != 0:
                    i, j = indices[a]
                    k, l = indices[b]
                    perms = supports[(i,j,k,l)]
                    new_mat = None
                    norm = None
                    for sigma in perms:
                        candidate = M - M[a,b] * compress(sigma, 2, n)
                        candidate_norm = np.linalg.norm(candidate, ord='fro')
                        if type(new_mat) == type(None):
                            new_mat = candidate
                            norm = np.linalg.norm(new_mat, ord='fro')
                        elif candidate_norm < norm:
                            new_mat = candidate
                            norm = candidate_norm
                    M = new_mat
                    for m in M:
                        print(list(m))
        counter += 1

        if counter == max_iter:
            # for m in M:
            #      print(list(m))
            print(c, d, e, f)
            print(':(')
            return 1
    return 0





find_coeff_4(1, 2, 2, 5, 5)
  
fail = []

# n = 5
# for i in range(1, n+1):
#     for j in range(1, n+1):
#         if i != j:
#             for k in range(1, n+1):
#                 for l in range(1, n+1):
#                     if k < l:
# #                        print(i,j,k,l)
#                         if find_coeff_4(i, j, k, l, n):
#                             fail.append((i, j, k, l))
print(fail)
print(len(fail))
n = 5
count = 0
k = 2
# for i in range(1, n+1):
#     for j in range(1, n+1):
#         if i != j:
#             for k in range(1, n+1):
#                 for l in range(1, n+1):
#                     if k != l:
#                         M = rep_w_ij_kl(i, j, k, l, n)
#                         for sigma in Permutation.group(n):
#                             if count_inc_seq(sigma, 2, n) > 0 and w_ij_kl(i, j, k, l)(sigma, n) == 1:
#                                 M = M -  compress(sigma, 2, n)
#                         if not (M == np.zeros(M.shape)).all():
#                             print(i, j, k, l)
#                             print(M)
#                             count += 1
                            

def list_to_freq_dict(l: list) -> dict:
    d = {}
    for i in l:
        if i in d.keys():
            d[i] += 1
        else:
            d[i] = 1
    return d

n = 5
i = 3
j = 5
k = 2
l = 1

M = rep_w_ij_kl(i, j, k, l, n)
# M = M - compress(Permutation(5, 3, 1, 2, 4), k, n)
# M = M - compress(Permutation(5, 3, 1, 2, 4), k, n)
# M = M - compress(Permutation(5, 3, 1, 2, 4), k, n)
# M = M - compress(Permutation(5, 3, 1, 2, 4), k, n)
# M = M + compress(Permutation(3, 4, 1, 2, 5), k, n)
# M = M + compress(Permutation(2, 3, 1, 5, 4), k, n)
# M = M + compress(Permutation(1, 3, 5, 2, 4), k, n)
# M = M - compress(Permutation(5, 3, 1, 2, 4), k, n)
# M = M + compress(Permutation(5, 2, 1, 3, 4), k, n)
# M = M + compress(Permutation(5, 1, 3, 2, 4), k, n)
# M = M + compress(Permutation(3, 4, 1, 2, 5), k, n)
# M = M - compress(Permutation(3, 2, 1, 4, 5), k, n)
# M = M + compress(Permutation(5, 2, 1, 3, 4), k, n)
# for m in M:
#     print(list(m))




# n = 6
# M = np.zeros((n**2 * (n - 1)**2, 78))
# col = 0
# perms = []

# n = 5
# for sigma in Permutation.group(n):
#     if count_inc_seq(sigma, k, n) > 0:
#         #print(compress(sigma, k, n).shape, n**2 * (n - 1)**2)
#         M[:, col] = np.reshape(compress(sigma, 2, n), (n**2 * (n - 1)**2, 1))[:, 0]
#         col += 1
#         perms.append(sigma)
# ls = np.linalg.inv((np.transpose(M) @ M)) @ (np.transpose(M))

# print(compress(Permutation(4, 1, 2, 3), 2, 4))

# for tup in [(1, 2, 4, 5), (1, 4, 4, 5), (2, 1, 4, 5), (2, 5, 3, 5), (3, 1, 2, 4), (3, 1, 2, 5), (3, 2, 2, 4), (3, 4, 1, 4), (3, 4, 1, 5), (3, 5, 2, 3), (4, 1, 3, 4), (4, 1, 3, 5), (4, 1, 4, 5), (4, 2, 3, 4), (4, 2, 3, 5), (4, 3, 1, 2), (4, 5, 1, 3), (5, 1, 2, 4), (5, 1, 2, 5), (5, 1, 3, 5), (5, 2, 3, 5)]:
#     print(tup)
#     v = ls @np.reshape(rep_w_ij_kl(*tup, n), (n**2 * (n - 1)**2, 1))
#     v = np.array([round(x, 5) for x in v[:,0]])
#     print({perms[i]: v[i] for i in range(len(v))})


# dist = {}
# for i in range(1, n+1):
#     for j in range(1, n+1):
#         if i != j:
#             for k in range(1, n+1):
#                 for l in range(1, n+1):
#                     if k < l:
# #                         print(i, j, k, l)
#                         v = ls @ np.reshape(rep_w_ij_kl(i, j, k, l, n), (n**2 * (n - 1)**2, 1))
#                         v = np.array([round(x, 5) for x in v[:,0]])
# #                         print(v)
#                         pairs = [(v[i], perms[i]) for i in range(len(v))]
#                         d = list_to_freq_dict(v)
#                         tup = 200*[0]
#                         for key, val in d.items():
#                             tup[100 + int(key)] = val
#                         tup = tuple(tup)
#                         if tup in dist.keys():
#                             dist[tup].append(i*1000 + j*100 + k*10 + l)
#                         else:
#                             dist[tup] = [i*1000 + j*100 + k*10 + l]
#                         print(dist[tup])
# print('-------------------------------------------------')                        
# for val in dist.values():
#     print(val)
#                         print(d)
#                         for a,b in pairs:
#                             if a != 0:
#                                 print(a,b, [b(i) for i in range(1, n+1)])
#                         print('\n')
# print(*perms)                        
                        # if not (np.reshape(M @ v, (n*(n-1), n*(n-1))) == rep_w_ij_kl(i, j, k, l, n)).all():
                        #     print(i, j, k, l)

#                         if d == {0: 72, 1: 6}:
#                             print(i, j, k, l)
#                             count += 1
# print(count)                        
# k = 4
# l = 5
# for i in range(1, n+1):
#     for j in range(1, n+1):
#         if i != j:
#             print(i, j, k, l)
#             v = ls @ np.reshape(rep_w_ij_kl(i, j, k, l, n), (n**2 * (n - 1)**2, 1))
#             v = np.array([round(x, 5) for x in v[:,0]])
#             print(v)
#             pairs = [(v[i], perms[i]) for i in range(len(v))]
#             d = list_to_freq_dict(v)
#             print(d)
#             for a,b in pairs:
#                 if a != 0:
#                     print(a,b, [b(i) for i in range(1, n+1)])
#             print('\n')


# n = 5
# M = np.zeros((20, 20))
# for sigma in Permutation.group(3):
#     tau = sigma * Permutation(1, 2, 5, 3, 4)
#     if tau in [Permutation(2, 1), Permutation(1, 3, 2), Permutation(3, 2, 1)]:
#         M = M + compress(tau, 2, n)
#     else:
#         M = M - compress(tau, 2, n)
# for m in M:
#     print(list(m))

# print()
# N = 2 * (compress(Permutation(1, 2, 5, 4, 3), 2, n) + compress(Permutation(1, 5, 2, 4, 3), 2, n) + compress(Permutation(5, 1, 2, 4, 3), 2, n))
# for n in N:
#     print(list(n))

# print()    
# for r in rep_w_ij_kl(4, 3, 4, 5, 5) - (M + N):
#     print(list(r))

n = 6
def find_full_support(n):
    for i in range(1, n+1):
        for j in range(1, n+1):
            if i != j:
                for k in range(1, n+1):
                    for l in range(1, n+1):
                        if k < l:
                            works = True
                            for sigma in Permutation.group(n):
                                if w_ij_kl(i, j, k, l)(sigma, n):
                                    if count_inc_seq(sigma, n, 2) == 0:
                                        works = False
                                        break
                            if works:
                                print(i, j, k, l)
# find_full_support(6)
                                    
                                
                            
# v = ls @ np.reshape(compress(Permutation(5, 4, 3, 1, 2), 2, n), (400, 1))
# one_hot = np.zeros((400, 1))
# one_hot[0][0] = 1
# v = ls @ one_hot
# v = np.array([round(600 * x, 5) for x in v[:,0]])

# print(v)
# n = 6
# i = 2
# j = 3
# k = 2
# l = 3
# print(i, j, k, l)
# v = ls @ np.reshape(rep_w_ij_kl(i, j, k, l, n), (n**2 * (n - 1)**2, 1))
# coeff= [round(x, 5) for x in v[:,0]]
# pairs = [(coeff[i], perms[i]) for i in range(len(v))]
# for a, b in pairs:
#     if a != 0:
#         print(a, ' ', [b(i) for i in range(1, n+1)])

# pos = np.array([max(i, 0) for i in coeff])
# neg = np.array([min(i, 0) for i in coeff])

# v_pos = np.reshape(M @ pos, (20, 20))
# v_neg = np.reshape(M @ neg, (20, 20))
# for row in v_pos:
#     print(list(row))
# for row in v_neg:
#     print(list(row))


# combo = []
# index = 0
# for sigma in Permutation.group(n):
#     if count_inc_seq(sigma, 2, n) > 0:
#         combo.append((str(sigma),  count_inc_seq(sigma, 2, n), coeff[index]))
#         index += 1
# print(combo)

# for a, b, c in combo:
#     if c == -2:
#         print(a, b)

# print(count)        

