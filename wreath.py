
def wreath(a: list) -> list:
    if len(a) == 1:
        return [a]
    else:
        left = wreath(a[0:int(len(a)/2)])
        right = wreath(a[int(len(a)/2):])
        result = []
        for x in left:
            for y in right:
                result.append(x + y)
                result.append(y + x)
        return result

print(wreath([1,2,3,4]))
    
