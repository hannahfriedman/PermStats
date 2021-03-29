from permutation import Permutation

def w_ij(sigma, j, i):
    """
    input: permutation
    output: 1 if sigma(j) = i
    """
    if sigma(j) == i:
        return 1
    else:
        return 0


def w_ij_kl(sigma, i, j, k, l):
    """
    input: permutation
    output: 1 if sigma(k) = i, sigma(l) = j 
    """
    if sigma(k) == i and sigma(l) == j:
        return 1
    else:
        return 0
