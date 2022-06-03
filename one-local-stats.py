import numpy as np
from itertools import combinations
from itertools import product
import matplotlib.pyplot as plt
from permutation import Permutation
import seaborn as sns
from D8_equivalence import D8_Equivalence

# positions = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]

# S3 = [np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]]), np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]), np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]]), np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])]

def generate_positions(n):
    return product(list(range(n)), list(range(n)))

positions = generate_positions(3)

def generate_reps(n):
    sn = Permutation.group(n)
    result = []
    for sigma in sn:
        mat = np.zeros((n,n))
        for col in range(1, n+1):
            mat[sigma(col) - 1, col - 1] = 1
        result.append(mat)
    return result

S3 = generate_reps(3)
S4 = generate_reps(4)

def powerset(iterable):
    s = list(iterable)
    pow_set = []
    for i in range(0, len(s) + 1):
        pow_set += combinations(s, i)
    return pow_set

def compute_dft_or(mats, indices):
    n = mats[0].shape[0]
    result = np.zeros((n, n))
    for mat in mats:
        count = 0
        for index in indices:
            if mat[index] == 1:
                count += 1
        result += count * mat
    return result

def compute_dft_and(mats, indices):
    n = mats[0].shape[0]
    result = np.zeros((n, n))
    for mat in mats:
        sets_of_indices = combinations(indices, 2)
        count = 0
        for pair in sets_of_indices:
            print(pair)
            print(mat[pair[0]], mat[pair[1]])
            if mat[pair[0]] == 1 and mat[pair[1]] == 1:
                count += 1
        result += count * mat
    return result

# def compute_dft_and_three(mats, indices):
#     n = mats[0].shape[0]
#     result = np.zeros((n, n))
#     sets_of_indices = combinations(indices, 3)
#     for mat in mats:
#         count = 0
#         for pair in sets_of_indices:
#             if mat[pair[0]] == 1 and mat[pair[1]] == 1 and mat[pair[2]] == 2:
#                 count += 1
#         result += count * mat
#     return result


def determine_real(vec):
    return np.isreal(np.round_(vec, 7)).all()

def determine_int(vec):
    return (np.round_(vec) == np.round_(vec, 7)).all()

def make_subplots(indices, nrows, ncols, n):
    for i in range(len(indices)):
        m = np.zeros((n, n))
        for index in indices[i]:
            m[index] = 1
        plt.subplot(nrows, ncols, i+1)
        plt.imshow(m, vmin=0, vmax=1)
        plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False, labeltop = False)
    plt.show()



## Plot Eigenvalues ##

# eigs = []
# for index in powerset(positions):
#     M = compute_dft_or(S3, index)
#     eig = np.linalg.eig(M)[0]
#     for val in eig:
#         eigs.append(val)
# eig_vec = np.array(eigs)
# x = eig_vec.real
# y = eig_vec.imag

# # plt.clf()
# sns.jointplot(x=x, y=y, kind='hex', color = 'purple')
# plt.show()

# heatmap, xedges, yedges = np.histogram2d(x, y, bins=75)
# extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]


# plt.imshow(heatmap.T, extent=extent, origin='lower', cmap="RdYlBu_r")
# plt.show()


# fig, ax = plt.subplots()
# plt.grid(color='gray', linestyle='-', linewidth=0.5)
# ax.set_xticks(list(range(-6, 20)))
# ax.scatter(eig_vec.real, eig_vec.imag, c=eig_vec.real, cmap="RdYlBu_r")
# plt.show()

## Partition according to real/integer ##

# real_int = []
# complex_int = []
# real_non_int = []
# complex_non_int = []
# zero = 0
# for index in powerset(positions):
#     M = compute_dft_and(S3, index)
#     if not M.any():
#         zero += 1
#     eig_vals = np.linalg.eig(M)[0]
#     if determine_real(eig_vals) and determine_int(eig_vals):
#         real_int.append(index)
#     elif determine_real(eig_vals):
#         real_non_int.append(index)
#     elif determine_int(eig_vals):
#         complex_int.append(index)
#     else:
#         complex_non_int.append(index)

# print(zero)
# print(len(real_int), len(complex_int), len(real_non_int), len(complex_non_int))

# make_subplots(complex_non_int, 3, 3, 3)

# make_subplots(real_int, 18, 19, 3)
# make_subplots(real_non_int, 9, 10, 3)
# make_subplots(complex_non_int, 10, 10, 3)

## Partition according to rank ##

ranks = [[], [], [], []]
eigs = []
for index in powerset(positions):
    M = compute_dft_and(S3, index)
    eig_vals = np.linalg.eig(M)[0]
    used = False
    for v in eigs:
        if (v == eig_vals).all():
            used = True
            break
    if not used:
        eigs.append(eig_vals)
    ranks[np.linalg.matrix_rank(M)].append(index)
    # if np.linalg.matrix_rank(M) == 3:
    #     print(np.linalg.eig(M)[0][0], len(index))
print(len(eigs))
# for l in ranks:
#     print(len(l))
for r in ranks:
    make_subplots(r, 16, 16, 3)
# make_subplots(ranks[1], 4, 4, 3)
# make_subplots(ranks[2], 15, 16, 3)
# make_subplots(ranks[3], 16, 17, 3)
# print([(0,0)] == [(0,0)] )

## Partition according to rank ##

# for index in powerset(positions):
#     mats = [index_to_matrix(index, 3)]
#     mats.append(np.rot90(mats[-1]))
#     mats.append(np.rot90(mats[-1]))
#     mats.append(np.rot90(mats[-1]))
#     mats.append(np.fliplr(mats[0]))
#     mats.append(np.flipud(mats[0]))
#     mats.append(np.transpose(mats[0]))
#     mats.append(np.transpose(mats[2]))
#     indices = [matrix_to_index(mat) for mat in mats]
#     eigs = np.linalg.eig(compute_dft_or(S3, index))[0]
#     real = determine_real(eigs)
#     integer = determine_int(eigs)
#     r = np.linalg.matrix_rank(compute_dft_or(S3, index))
#     print('---------')
#     show = False
#     for i in range(1, 8):
#         if np.linalg.matrix_rank(compute_dft_or(S3, indices[i])) != r:
#             print(i)
#             show = True
#     if show:
#         make_subplots(indices, 2, 4, 3)            
#         if determine_real(np.linalg.eig(compute_dft_or(S3, indices[i]))[0]) != real or determine_int(np.linalg.eig(compute_dft_or(S3, indices[i]))[0]) != integer:
#             print(i)
#             show = True


## Plot eigenvalues by rotation ##
count = 1
for equiv_class in D8_Equivalence.powerset_equivalence(positions):
    print(equiv_class.equivalent)
    eigs = []
    for indices in equiv_class.equivalent:
        for eig in np.linalg.eig(compute_dft_or(S3, indices))[0]:
            eigs.append(eig)
#     plt.subplot(10, 11, count)
    sns.jointplot(x=np.array(eigs).real, y=np.array(eigs).imag, kind='hex', color = 'purple')
    count += 1
    plt.show()
    


       


# for i in range(len(real_int)):
#     m = np.zeros((3, 3))
#     for index in real_int[i]:
#         m[index] = 1
#     plt.subplot(18, 19, i+1)
#     plt.imshow(m, vmin=0, vmax=1)
#     plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False, labeltop = False)
# plt.show()

# for i in range(len(real_non_int)):
#     m = np.zeros((3, 3))
#     for index in real_non_int[i]:
#         m[index] = 1
#     plt.subplot(9, 10, i+1)
#     plt.imshow(m, vmin=0, vmax=1)
#     plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False, labeltop = False)
# plt.show()

# for i in range(len(complex_non_int)):
#     m = np.zeros((3, 3))
#     for index in complex_non_int[i]:
#         m[index] = 1
#     plt.subplot(10, 10, i+1)
#     plt.imshow(m, vmin=0, vmax=1)
#     plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False, labeltop = False)
# plt.show()

# print(real_int[0], real_int[-1])

        
    
