import unittest

from BaseTree import BaseTernaryTree, QubitNum, BranchNum, NodeContacts, LostNum
from TernaryTree import TernaryTree, pauli_weight, prod_pauli_strings


from TreeVizualization import draw


def bonsai():
    # tt = TernaryTree(nodes=[(0,),
    #                         (1,2,3),
    #                         (None,None,4, None, None, 5, None, None, 6),
    #                         (7, None,8, 9, None,10, 11,None,12),
    #                         (None,None,13,None,None,14,None,None,15,None,None,16,None,None,17,None,None,18),
    #                         (19,None,20,21,None,22,23,None,24,25,None,26,27,None,28,29,None,30) ,
    #                         (None,None,31,None,None,32,None,None,33,None,None,34,None,None,35,None,None,36,None,None,37,None,None,38,None,None,
    #                          39,None,None,40,None,None,41,None,None,42),
    #                         (43,None,44,45,None,46,47,None,48,49,None,50,51,None,52,53,None,54,55,None,56,57,None,58,59,None,60,61,None,62,63,None,64,65,None,66,
    #                          67,None,68,69,None,70,71,None,72,73,None,74,75,None,76,77,None,78)])
    tt = TernaryTree(nodes=[(0,),
                            (1, None, 2),
                            (3, 4, 5, 6, 7, 8),
                            (None, None, 13, None, None, 14, None, None, 15, None, None, 16, None, None, 17, None, None,
                             18),
                            (19, None, 20, 21, None, 22, 23, None, 24, 25, None, 26, 27, None, 28, 29, None, 30),
                            (None, None, 31, None, None, 32, None, None, 33, None, None, 34, None, None, 35, None, None,
                             36, None, None, 37, None, None, 38, None, None,
                             39, None, None, 40, None, None, 41, None, None, 42),
                            (43, None, 44, 45, None, 46, 47, None, 48, 49, None, 50, 51, None, 52, 53, None, 54, 55,
                             None, 56, 57, None, 58, 59, None, 60, 61, None, 62, 63, None, 64, 65, None, 66,
                             67, None, 68, 69, None, 70, 71, None, 72, 73, None, 74, 75, None, 76, 77, None, 78)
                            ])
    tt.update_branchnum()

    return tt


def st_enumeration1(nmodes: int = 36):
    """
    Standard numeration for alpha_beta_tree which is effective on my opinion
    """
    num_list = []
    j = 0
    if nmodes % 2 == 0 and nmodes > 2:
        for i in range(1, nmodes - 1):
            if j <= 1:
                j += 1
                num_list = num_list + [i, i + 2]
            else:
                if j == 2:
                    j = 3
                else:
                    j = 0
        if nmodes //2 % 2 != 0:
            num_list = num_list + [nmodes -1, nmodes]
        for i in range(1, nmodes - 1):
            if j <= 1:
                j += 1
                num_list = num_list + [ i + nmodes, i + 2 + nmodes]
            else:
                if j == 2:
                    j = 3
                else:
                    j = 0
        if nmodes //2 % 2 != 0:
            num_list = num_list + [nmodes*2 - 1, 2*nmodes ]
    return num_list


def st_enumeration2(nmodes: int = 36):
    """
    Standard numeration for alpha_beta_tree which is effective on my opinion
    """
    num_list = []

    if nmodes % 2 == 0 and nmodes > 2:
        j = 1
        for i in range(1, 18):
            num_list = num_list + [i*j + (nmodes - i + 1) * (1-j)]
            j = (j + 1) % 2
        j = 1
        for i in range(1, 19):
            num_list = num_list + [(i + 1)*j + (nmodes - i + 2) * (1-j)]
            j = (j + 1) % 2
        j = 1
        for i in range(37, 55):
            num_list = num_list + [i * j + (3*nmodes - i + 1) * (1 - j)]
            j = (j + 1) % 2
        j = 1
        for i in range(37, 55):
            num_list = num_list + [(i + 1) * j + (3*nmodes - i + 2) * (1 - j)]
            j = (j + 1) % 2

    return num_list


def st_enumeration3(nmodes: int = 30):
    """
    Standard numeration for alpha_beta_tree which is effective on my opinion
    """
    # nmodes = 36
    k = nmodes // 4 + 1
    j = nmodes // 2 % 2
    num_list = []

    if nmodes % 2 == 0 and nmodes > 2:
        for i in range(1, k):
            num_list = num_list + [2 * i - 1] + [nmodes - 2 * i + 1]
        if j == 1:
            num_list = num_list + [nmodes//2]

        for i in range(1, k):
            num_list = num_list + [2 * i ] + [nmodes - 2 * i + 2]
        if j == 1:
            num_list = num_list + [nmodes//2 + 1]

        for i in range(1, k):
            num_list = num_list + [nmodes + 2 * i - 1] + [2*nmodes - 2 * i + 1]
        if j == 1:
            num_list = num_list + [nmodes//2 + nmodes]

        for i in range(1, k):
            num_list = num_list + [nmodes + 2 * i] + [2*nmodes - 2 * i + 2]
        if j == 1:
            num_list = num_list + [nmodes//2 + nmodes + 1]

    return num_list


class TestNode(unittest.TestCase):
    def setUp(self):
        self.q = QubitNum()
        self.b = BranchNum()

    def test_equations(self):
        num = 1
        childs = None
        q = NodeContacts(1)
        self.assertEqual(isinstance(q.parent, QubitNum), True)
        self.assertEqual(q.parent.num, num,)
        for child in q.childs:
            self.assertEqual(isinstance(child, LostNum), True)
        num = -1
        childs = [1, 2]
        q = NodeContacts(num, childs)
        self.assertEqual(isinstance(q.parent, QubitNum), True)
        # self.assertEqual(isinstance(q.parent, QubitNum), True)

        for i, child in enumerate(childs):
            self.assertEqual(q.childs[i], QubitNum(child))
        self.assertEqual(isinstance(q.childs[2], LostNum), True)

        q[2] = 4
        self.assertEqual(q.childs[2], QubitNum(4))
        # self.assertEqual(0 == 1, True)


class TestNodeNum(unittest.TestCase):

    def test_equations(self):
        q = QubitNum(1)
        b = BranchNum(1)
        q1 = QubitNum(1)
        q2 = QubitNum(10)
        b1 = BranchNum(1)
        b2 = BranchNum(10)
        # self.assertEqual(q == b, False, "objects of different classes  are the same")
        # self.assertEqual(b == q, False, "objects of different classes  are the same")
        self.assertEqual(q1 == q2, False, "Different QubitNum objects equal")
        self.assertEqual(b1 == b2, False, "Different BranchNum objects equal")
        self.assertEqual(None == 0, False)
        self.assertEqual(q1 < q2, True)
        self.assertEqual(q1 >  q2, False)
        self.assertEqual(q1 < q1, False)
        # self.assertEqual(q1 < b1, False)
        self.assertEqual(b1 < b2, True)
        self.assertEqual(b1 > b2, False)
        self.assertEqual(b1 < b1, False)
        k = {}
        a = QubitNum(1)
        b = QubitNum(1)
        k[a] = 1
        k[b] = 2
        self.assertEqual(k[a] == 1, False)
        self.assertEqual(k[b] == 2, True)


class TestTernaryTree(unittest.TestCase):

    def test_equations(self):
        tt = TernaryTree(1)
        tt = TernaryTree(5)
        # print(tt)
        for branch in tt.branches():
            # print(branch)
            pass
        tt.update_branchnum(list(range(10, 0, -1)))
        # print(tt)
        s = tt.branch_transposition(1,2,0,1)
        # print("s = ", s)
        # print(tt)
        tt = TernaryTree(10)
        # print(tt.nodes)
        tt.delete_node(5)
        # print(tt.nodes)сщтвф
        tt.add_node(5, 4, 0)
        # print(tt.nodes)
        tt.update_branchnum()
        # print(tt.nodes)
        tt[5][2] = LostNum()
        # print(tt.nodes)
        tt = TernaryTree(nodes = [(0,),(1,2,3),(4,5,None,6,None, 7,8,9)])
        # print(tt.nodes)

    def test_bonsai(self):
        pass
        # tt = bonsai()
        # # print(tt.nodes)
        # # print(tt.num_nodes)
        # tt.update_branchnum()
        # # print(tt.nodes)
        # tt[36][2] = LostNum()
        #
        # # print(tt)

    def test_draw(self):
        # tt = TernaryTree(4)
        # tt = bonsai()

        # draw(tt)
        pass

    def test_prod(self):
        import numpy as np
        nmodes = 40
        # tt = bonsai()
        qubits = list(range(0,nmodes))
        t = []
        l = 0
        while 3**l + (3**l - 1)//2 < nmodes:

            t = t + [(qubits[(3**l - 1)//2:3**l + (3**l - 1)//2])]
            l += 1

        t = t + [(qubits[(3 ** l - 1) // 2:nmodes])]
        # print(t)
        tt = TernaryTree(nodes = t)
        # tt.delete_node(66)
        # tt.delete_nodes([34,35,36])
        tt.update_branchnum()

        # print(tt.branches())
        tt[39][2] = LostNum()
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
            occupied = 9
            par = (i + j) % 2 == 1
            # beta = (i>36 and j>36 and (i <= 2*occupied + 36) and (j > 2*occupied + 36))
            # alpha = (i <= 2*occupied) and (j > 2*occupied) and i<=36 and j<=36
            # return par and (alpha or beta)
            return par

        def valid_4prod(i,j,k,l):
            occupied = 9
            par = (i+j+k+l)%2 == 0
            # alphaij = (i <= 2 * occupied) and (j > 2 * occupied) and i <= 36 and j <= 36
            # betaij = (i > 36 and j > 36 and (i <= 2 * occupied + 36) and (j > 2 * occupied + 36))
            # alphakl = (k <= 2 * occupied) and (l > 2 * occupied) and k <= 36 and l <= 36
            # betakl = (k > 36 and l > 36 and (k <= 2 * occupied + 36) and (l > 2 * occupied + 36))
            # return par and (alphaij and alphakl or alphaij and betakl or betaij and betakl)
            return par

        draw(tt,k = 1)
        s = 0
        occupied = 9

        # print(pauli_weight(tt, [25, 7]))
        n = 0
        for i in range(1, 2*nmodes + 1):
            if i <= nmodes:
                m1 = nmodes + 1
            else:
                m1 = 2*nmodes + 1
            for j in range(i + 1, m1):
                if valid_2prod(i,j):
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


        print(s)
        print(n)
        print(s/n)
if __name__ == "__main__":
    unittest.main()