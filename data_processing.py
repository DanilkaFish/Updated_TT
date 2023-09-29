from __init__ import *


def main():
    nmodes = 40
    tt = TernaryTree()
    tt.build_full_tt(nmodes)

    # tt.update_branchnum()

    # print(tt.branches())
    # tt[39][2] = LostNum()
    draw(tt)
    tt.to0vac()

    # draw(tt)
    import numpy as np
    # enum_list = [int(i) for i in np.random.permutation(list(range(1, 61)))]
    # enum_list = st_enumeration3(nmodes)
    # enum_list = [1,5,2,4,6,3,7,9,11,12,8,10]
    # print(enum_list)
    # enum_list = list(range(1, nmodes*2 + 1))
    # enum_list = [1,2,3,4,5,6]
    # # enum_list = list(range(1,26)) + list(range(27,51)) + [26] + list(range(51,73))
    # tt.update_branchnum(enum_list)

    def valid_2prod(i, j):
        par = (i + j) % 2 == 1
        return par

    def valid_4prod(i, j, k, l):
        par = (i+j+k+l)%2 == 0
        return par

    draw(tt, k = 1)

    s = 0
    n = 0

    for i in range(1, 2*nmodes + 1):
        if i <= nmodes:
            m1 = nmodes + 1
        else:
            m1 = 2*nmodes + 1
        for j in range(i + 1, m1):
            if valid_2prod(i, j):
                r = pauli_weight(tt, [i, j])
                if r > 0:
                    s += r
                    n += 1

    for i in range(1, 2*nmodes + 1):
        if i <= nmodes :
            m1 = nmodes + 1
        else:
            m1 = 2*nmodes + 1
        for j in range(i, m1):
            for k in range(j, 2*nmodes + 1):
                if i <= nmodes:
                    m2 = nmodes + 1
                else:
                    m2 = 2*nmodes + 1
                for l in range(k, m2):
                    if valid_4prod(i, j, k, l):
                        r = pauli_weight(tt, [i, j, k, l])
                        if r >0:
                            s += r
                            n += 1


if __name__ == "__main__":
    main()