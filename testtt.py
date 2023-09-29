import unittest

# from BaseTree import BaseTernaryTree, QubitNum, BranchNum, NodeContacts, LostNum
# from TernaryTree import TernaryTree, pauli_weight, prod_pauli_strings
from __init__ import *

from TreeVizualization import draw


def bonsai():
    flag = 0
    if flag:
        tt = TernaryTree(nodes=[(0,),
                                (1,2,3),
                                (None,None,4, None, None, 5, None, None, 6),
                                (7, None,8, 9, None,10, 11,None,12),
                                (None,None,13,None,None,14,None,None,15,None,None,16,None,None,17,None,None,18),
                                (19,None,20,21,None,22,23,None,24,25,None,26,27,None,28,29,None,30) ,
                                (None,None,31,None,None,32,None,None,33,None,None,34,None,None,35,None,None,36,None,None,37,None,None,38,None,None,
                                 39,None,None,40,None,None,41,None,None,42),
                                (43,None,44,45,None,46,47,None,48,49,None,50,51,None,52,53,None,54,55,None,56,57,None,58,59,None,60,61,None,62,63,None,64,65,None,66,
                                 67,None,68,69,None,70,71,None,72,73,None,74,75,None,76,77,None,78)])
    else:
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
        self.q = QubitNum(1)
        self.b = BranchNum(1)
        self.l = LostNum()

    def test_equations(self):
        num = 1
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

    def test_last(self):
        self.assertTrue(not self.l)
        self.assertTrue(self.b)
        self.assertTrue(self.q)
        self.assertTrue(self.b.is_last)
        self.assertTrue(not self.q.is_last)
        self.assertTrue(self.l.is_last)


class TestNodeNum(unittest.TestCase):

    def test_equations(self):
        q = QubitNum(1)
        b = BranchNum(1)
        q1 = QubitNum(1)
        q2 = QubitNum(10)
        b1 = BranchNum(1)
        b2 = BranchNum(10)

        self.assertEqual(q1 == q2, False, "Different QubitNum objects equal")
        self.assertEqual(b1 == b2, False, "Different BranchNum objects equal")
        self.assertEqual(None == 0, False)
        self.assertEqual(q1 < q2, True)
        self.assertEqual(q1 >  q2, False)
        self.assertEqual(q1 < q1, False)
        self.assertEqual(q1 < b1, False)
        self.assertEqual(b1 < b2, True)
        self.assertEqual(b1 > b2, False)
        self.assertEqual(b1 < b1, False)

    def test_hash(self):
        k = {}
        a = QubitNum(1)
        b = QubitNum(1)
        c = BranchNum(1)
        k[c] = 3
        k[b] = 2
        k[a] = 1
        self.assertEqual(k[a], 1)
        self.assertNotEqual(k[b], 2)
        self.assertEqual(k[c], 3)


class TestTernaryTree(unittest.TestCase):

    def test_initialization(self):
        tt = TernaryTree()
        tt = TernaryTree(0)
        tt = TernaryTree(5)
        tt.update_branchnum(list(range(10, 0, -1)))
        tt.branch_transposition(1,2,0,1)

        tt.build_full_tt(10)
        tt.build_full_tt(0)
        tt.build_full_tt(4)
        tt = TernaryTree(10)
        tt.delete_node(5)

        tt.add_node(5, 4, 0)
        tt.delete_node([5,6,7,8,9])
        tt.update_branchnum()
        draw(tt)
        tt.update_branchnum(enum_list = range(1, 2*tt.num_nodes + 2))
        self.assertRaises(TreeStructureError, lambda: tt.update_branchnum(enum_list = range(1, 2)))
        tt.change_num(4, 2)
        tt = TernaryTree(nodes = [(0,),(1,2,3),(4,5,None,6,None, 7,8,9)])

    def test_draw(self):
        tt = TernaryTree(4)
        # tt = bonsai()
        # draw(tt)

    def test_to0vac(self):
        tt = TernaryTree(10)
        tt.build_full_tt(5)
        tt.change_num(4,2,11)
        self.assertRaises(TreeStructureError, tt.to0vac)
        tt.change_num(4,2)
        tt.change_num(4,1)

        self.assertRaises(IndexError, tt.to0vac)
        # tt[4][2] = LostNum()
        # tt.to0vac()


if __name__ == "__main__":
    unittest.main()