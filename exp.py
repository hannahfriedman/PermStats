from permutation import Permutation
from misc import *

class Tup(object):
    n = 5
    def __init__(self, a: int, b: int):
        self.data = (a, b)
        if a > b:
            self.bt = a - b -1
            self.lt = b - 1
            self.gt = self.n - a
            self.df = [(a-b-1, self.n - a), (b-1, a - b - 1)]
        if b > a:
            self.bt = b - a - 1
            self.lt = a - 1
            self.gt = self.n - b
            self.df = [(a - 1, b - a - 1), (b - a - 1, self.n - b)]

    def __repr__(self):
        return str(self.data)
    
    @classmethod
    def gen_all_inc(cls):
        result = []
        for i in range(1, cls.n+1):
            for j in range(i+1, cls.n+1):
                result.append(Tup(i, j))
        return result

    @classmethod
    def gen_all_dec(cls):
        result = []
        for i in range(1, cls.n+1):
            for j in range(1, i):
                result.append(Tup(i, j))
        return result

    @classmethod
    def print_pairs2_fix3(cls):
        all_inc = Tup.gen_all_inc()
        all_dec = Tup.gen_all_dec()
        total = 0
        for ind in all_inc:
            for img in all_dec:
                contr = (ind.bt + ind.lt) * (img.bt + img.gt) + (ind.bt + ind.gt) * (img.bt + img.lt) + ind.lt * img.gt + ind.gt * img.lt
                contr1 = sum(ind.df[0]) * sum(img.df[0]) + sum(ind.df[1]) * sum(img.df[0]) + ind.df[0][0] * img.df[0][1] + ind.df[1][1] * img.df[1][0]
                print(f'{ind} -> {img}: {contr}')
                print(contr1)
                total += contr
                if contr != contr1:
                    print("ind: ", ind.lt, ind.bt, ind.gt, ind.df)
                    print("img: ", img.gt, img.lt, img.lt, ind.df)
        print("Total: " + str(total/2))

    @classmethod
    def print_pairs3_fix3(cls):
        all_inc = Tup.gen_all_inc()
        all_dec = Tup.gen_all_dec()
        total = 0
        for ind in all_inc:
            for img in all_dec:
                contr = img.bt * ind.bt # img.df[0][0] * ind.df[1][0] # * (ind.df[1][0] + ind.df[0][1])
                print(f'{ind} -> {img}: {contr}')
                total += contr
        print("Total: " + str(total))


    @classmethod
    def print_pairs3_fix4(cls):
        all_inc = Tup.gen_all_inc()
        all_dec = Tup.gen_all_dec()
        total = 0
        for ind in all_inc:
            for img in all_dec:
                contr = sum(ind.df[0]) * sum(img.df[0]) + sum(ind.df[1]) * sum(img.df[0]) + ind.df[0][0] * img.df[0][1] + ind.df[1][1] * img.df[1][0]
                #                print(f'{ind} -> {img}: {contr}')
                total += contr
                print("Total: " + str(total/2))


        
if __name__ == "__main__":
    Tup.print_pairs2_fix3()
'''
S4: 64 = 8^2
S5: 400 = 20^2
S6: 1600 = 40^2
S7: 4900 = 70^2
S8: 12544 = 112^2
S9: 168^2
Sn: (2 * (n) choose 3)^2

S4: 80 = 64 + 16 = 8^2 + 4^2
S5: 500 = 400 + 100 = 20^2 + 10^2
S6: 2000 = 1600 + 400 = 40^2 + 20 ^2
S7: 6125 = 4900 + 1225 = 70^2 + 35^2
S8: 15680 = 12544 + 3136 = 112^2 + 56^2
S9: 35280 = 168^2 + 84^2
S10: 72000 = 240^2 + 120^2
S11: 136125 = 330^2  + 165^2  
S12: 242000 = 440^2 + 220^2
S13: 408980 = 572^2  + 286^2
Sn: 5*(n choose 3)^2  
'''
