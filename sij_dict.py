from permutation import Permutation

class Sij_dict(object):
    def __init__(sigma: Permutation, k: int, n: int) -> None:
        self.data = {}

    class Sij(object):
        def __init__(self, ind: tuple, img: tuple, parents=[]):
            self.ind = ind
            self.img = img
            self.parents = parents
    
        def __getitem__(self, index: int) -> int:
            try:
                return self.img[self.ind.index(index)]
            except Exception:
                raise Exception(f"{ind}'s mapping is not specified in {self}.")
    
        def __repr__(self):
            return str(self.ind) + " -> " + str(self.img)
  
        def __eq__(self, other) -> bool:
            return self.ind == other.ind and self.img == other.img
    
        def __len__(self):
            return len(self.ind)
    
        def __hash__(self):
            return hash(repr(self))
