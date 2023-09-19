import unittest

from BaseTree import BaseTernaryTree, QubitNum, BranchNum, NodeContacts, LostNum
from TernaryTree import TernaryTree, pauli_weight, prod_pauli_strings


from TreeVizualization import draw


def bonsai():
    tt = TernaryTree(nodes=[(0,),
                            (1,2,3),
                            (None,None,4, None, None, 5, None, None, 6),
                            (7, None,8, 9, None,10, 11,None,12),
                            (None,None,13,None,None,14,None,None,15,None,None,16,None,None,17,None,None,18),
                            (19,None,20,21,None,22,23,None,24,25,None,26,27,None,28,29,None,30) ,
                            (None,None,None,None,None,None,None,None,None,None,None,31,None,None,None,None,None,None,None,None,None,None,None,32,None,None,
                             None,None,None,None,None,None,None,None,None,33,None,None),
                            (None,None,34,None,None,35,None,None,36)])
    tt.update_branchnum()

    return tt


def st_enumeration1(nmodes: int = 36):
    """
    Standard numeration for alpha_beta_tree which is effective on my opinion
    """
    num_list = []
    if nmodes % 2 == 0 and nmodes > 2:
        for i in range(1, 1 + nmodes // 2):
            num_list = num_list + [i, nmodes - i + 2 *(1 -  i % 2)]
        for i in range(nmodes + nmodes, nmodes + nmodes // 2 , -1):
            num_list = num_list + [i, - i + nmodes * 3 + 2 *(1 -  i % 2)]
    return num_list

def st_enumeration2(nmodes: int = 36):
    """
    Standard numeration for alpha_beta_tree which is effective on my opinion
    """
    num_list = []

    if nmodes % 2 == 0 and nmodes > 2:
        j = 1
        for i in range(1, 19):
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

        tt = bonsai()
        tt.delete_node(36)
        tt[33][2] = LostNum()
        import numpy as np
        # enum_list = [int(i) for i in np.random.permutation(list(range(1, 73)))]
        # enum_list = st_enumeration2()
        # print(enum_list)
        # enum_list = list(range(1, 73))
        enum_list = list(range(1,26)) + list(range(27,51)) + [26] + list(range(51,73))
        tt.update_branchnum(enum_list)

        def valid_2prod(i, j):
            occupied = 9
            par = (i + j) % 2 == 0
            beta = (i>36 and j>36 and (i <= 2*occupied + 36) and (j > 2*occupied + 36))
            alpha = (i <= 2*occupied) and (j > 2*occupied) and i<=36 and j<=36
            return par and (alpha or beta)

        def valid_4prod(i,j,k,l):
            occupied = 9
            par = (i+j+k+l)%2 == 0
            alphaij = (i <= 2 * occupied) and (j > 2 * occupied) and i <= 36 and j <= 36
            betaij = (i > 36 and j > 36 and (i <= 2 * occupied + 36) and (j > 2 * occupied + 36))
            alphakl = (k <= 2 * occupied) and (l > 2 * occupied) and k <= 36 and l <= 36
            betakl = (k > 36 and l > 36 and (k <= 2 * occupied + 36) and (l > 2 * occupied + 36))
            return par and(alphaij and alphakl or alphaij and betakl or betaij and betakl)
        draw(tt)
        s = 0
        occupied = 9

        print(pauli_weight(tt, [25, 7]))
        n = 0
        for i in range(1, 73):
            for j in range(1, 73):
                if valid_2prod(i,j):
                    s += pauli_weight(tt, [i, j])
                    n += 1
        for i in range(1,73):
            for j in range(1, 73):
                for k in range(1, 73):
                    for l in range(1, 73):
                        if valid_4prod(i, j, k, l):
                            n += 1
                            s += pauli_weight(tt, [i, j, k, l])
        print(s)
        print(n)

if __name__ == "__main__":
    unittest.main()