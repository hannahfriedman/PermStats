def statistic(pi):
    if pi.cycle_type() == [1,1,1]:
        return 12
    if pi.cycle_type() == [2,1]:
        return 6
    if pi.cycle_type() == [3]:
        return 3
    if pi.cycle_type() == [1,1,1,1]:
        return 48
    if pi.cycle_type() == [2,1,1]:
        return 32
    if pi.cycle_type() == [3,1]:
        return 24
    if pi.cycle_type() == [2,2]:
        return 16
    if pi.cycle_type() == [4]:
        return 16
    if pi.cycle_type() == [1,1,1,1,1]:
        return 240
    if pi.cycle_type() == [2,1,1,1]:
        return 180
    if pi.cycle_type() == [3,1,1]:
        return 150
    if pi.cycle_type() == [2,2,1]:
        return 120
    if pi.cycle_type() == [4,1]:
        return 120
    if pi.cycle_type() == [3,2]:
        return 90
    if pi.cycle_type() == [5]:
        return 90



