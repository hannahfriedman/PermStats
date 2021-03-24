from permutation import Permutation
def wolfram_syntax(mat):
    string = "{"
    for row in range(len(mat)):
        string+="{"
        for col in range(len(mat[0])-1):
            string += str(mat[row, col])
            string += ","
        string += str(mat[row, len(mat[0])-1])
        if row == len(mat) - 1:
            string += "}"
        else:
            string += "},"
    string+= "}"
    print(string)

def matlab_syntax(mat):
    string = "["
    for row in range(len(mat)):
        for col in range(len(mat)):
            string += " "
            string += str(mat[row, col])
        if row == len(mat) - 1:
            string += "]"
        else:
            string += ";"
    print(string)

