from sage.combinat.permutation import Permutation


def oij(I: list, J: list):
    def f(sigma: Permutation, n: int):
        for i in range(len(I)):
            if sigma(I[i]) != J[i]:
                return 0
        return 1
    return f
