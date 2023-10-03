import unittest

# from BaseTree import BaseTernaryTree, QubitNum, BranchNum, NodeContacts, LostNum
# from TernaryTree import TernaryTree, pauli_weight, prod_pauli_strings
from __init__ import *

# from TreeVizualization import draw


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
        tt.build_bk_tree(10)
        tt.build_bk_tree(0)
        tt.build_bk_tree(4)
        tt.build_alphabeta_tree(4)
        tt.build_alphabeta_tree(8)
        tt.build_alphabeta_tree(14)
        # print(tt.get_majorana(7))
        # print(tt.get_majorana(28))
        self.assertRaises(ValueError, lambda: tt.get_majorana(0))
        tt = TernaryTree(10)
        tt.delete_node(5)

        tt.add_node(5, 4, 0)
        tt.delete_node([5,6,7,8,9])
        tt.update_branchnum()
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