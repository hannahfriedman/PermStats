def matrix_rep(n):
    partitions = generate_partitions(n)
    tableaux_by_shape = [tableaux_shape(n, partition) for partition in partitions]
    rho = {}
    for i in range(1,n):
        representation = []
        for shape in tableaux_by_shape:
            sort_tableaux(n, shape)
            print(shape)    
        rho["(" + str(i) + "," + str(i+1) + ")"] = representation
    return rho


def generate_partitions(n):
    ans = []
    used = []
    if n == 1:
        return [[1]]
    elif n == 0:
        return [[]]
    for x in range(1, n):
        ans+=[[x] + part for part in generate_partitions(n-x)]
    return remove_dubs(ans) + [[n]]

def remove_dubs(partition):
    for part in partition:
        part.sort()
        part.reverse()
    return rem_dubs_helper(partition)

def rem_dubs_helper(partition):
    if len(partition) == 1:
        return partition
    for i in range(1, len(partition)):
        if partition[0] == partition[i]:
            return rem_dubs_helper(partition[1:])
    return [partition[0]] + rem_dubs_helper(partition[1:])


def generate_tableaux(n):
    if n == 1:
        return [[[1]]]
    else:
        prev_tab = generate_tableaux(n-1)
        ans = [tab + [[n]] for tab in prev_tab]
        for tab in prev_tab:
            for i in range(0, len(tab)):
                if len(tab) == 1:
                    ans += [[tab[i] + [n]]]
                elif i == 0 or len(tab[i-1]) > len(tab[i]):
                    ans += [tab[0:i] + [tab[i] + [n]] + tab[i+1:]]
        return ans

def tableaux_shape(n, partition):
    ans = []
    matching = True
    tableaux = generate_tableaux(n)
    for tab in tableaux:
        if len(tab) == len(partition):
            for row_index in range(len(tab)):
                if len(tab[row_index]) != partition[row_index]:
                    matching = False
            if matching:
                ans += [tab]
            else:
                matching = True
    return ans
                
def print_tableaux(tableau):
    for row in tableau:
        print(row)
    return

def compare_tableaux(n, tableau1, tableau2):
    for row_index in range(len(tableau1)):
        if n in tableau1[row_index] and n in tableau2[row_index]:
            return compare_tableaux(n-1, tableau1, tableau2)
        elif n == tableau1[row_index][-1]:
            return True
        elif n == tableau2[row_index][-1]:
            return False

def sort_tableaux(n, tableaux):
    switched = True
    while switched:
        switched = False
        for i in range(len(tableaux) - 1):
            if compare_tableaux(n, tableaux[i+1], tableaux[i]):
                tableaux[i], tableaux[i+1] = tableaux[i+1], tableaux[i]
                switched = True
    return


def 

