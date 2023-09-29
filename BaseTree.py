from __future__ import annotations

import warnings

from Nodes import BranchNum, QubitNum, LostNum, NodeContacts, ROOT


gate_name_to_number = {'X': 0, 'Y': 1, 'Z': 2}
gate_number_to_name = {0: 'X', 1: 'Y', 2: 'Z'}
dict_prod = {"II": "I", "XX": "I", "YY": "I", "ZZ": "I",
             "XY": 'Z', "YX": 'Z', "XZ": 'Y', "ZX": 'Y', "YZ": 'X', "ZY": 'X',
             "IX": "X", "XI": "X", "YI": "Y", "IY": "Y", "IZ": "Z", "ZI": "Z"}


class TreeStructureError(Exception):
    def __init__(self, text):
        self.txt = text


def st_enumeration(nmodes):
    """
    Standard numeration for alpha_beta_tree which is effective on my opinion
    """
    num_list = []
    if nmodes % 2 == 0 and nmodes > 2:
        for i in range(nmodes//2):
            num_list.append([i, nmodes - i - 1])
        for i in range(nmodes + nmodes - 1, nmodes + nmodes//2 - 1, -1):
            num_list.append([i, - i + nmodes*3 - 1])
    else:
        num_list = [[2*i, 2*i + 1] for i in range(nmodes)]
    return num_list


class BaseTernaryTree:
    """ This is the class for representation ternary tree and its manipulations"""
    def __init__(
            self,
            n_qubits: int = 0,
            nodes: list[tuple[int]] = None,
            root: QubitNum | int = QubitNum(0)
    ):
        """
        self.n_qubits = int -- number of qubits or number of tree nodes
        self.nodes = {parent_number : NodeContacts} -- tree's structure
        self.num_particles = tuple(int,int) -- number alpha and beta electrons respectively (Need only for branches' numeration)
        self.nmodes = int number -- number fermionic modes (Need only for branches' numeration)
        self.enum_list = [int] -- branches numeration
        """
        self.root = QubitNum(root)
        self._nodes = {}
        if nodes is not None:
            self.nodes = nodes
        else:
            self.build_st_tree(n_qubits)

    @property
    def nodes(self):
        return self._nodes

    def parent(self, node: QubitNum | int):
        return self[node].parent

    def change_num(self, node: int | QubitNum, edge: int, num: int | BranchNum = None):
        node = QubitNum(node)

        if num is None:
            num = LostNum()
        else:
            num = BranchNum(num)
        try:
            self[node][edge] = num
        except KeyError:
            warnings.warn("Here is no node ", node)
        try:
            self.check_branch_numeration()
        except TreeStructureError:
            self.update_branchnum(renumerate=False)

    def is_root(self, node: QubitNum | int):
        if self[node].parent == ROOT:
            return True
        else:
            return False

    @property
    def num_nodes(self):
        return len(self.nodes)

    def height(self):
        """
        Return tree's height
        """
        branches = self.branches()
        height = max([len(branch[1]) for branch in branches])
        return height

    @nodes.setter
    def nodes(self, nodes):
        self._nodes = {}
        last_level = [None]
        max_level_nodes = 1
        for index, level in enumerate(nodes):
            minim = min(len(level), max_level_nodes)
            number_level_nodes = 0
            for i in range(minim):
                if isinstance(level[i], int):
                    number_level_nodes += 1
                    self.add_node(node=level[i], its_parent=last_level[i//3], pos=i % 3)

            last_level = [i for i in level if not (i is None)]
            max_level_nodes = number_level_nodes*3

    def build_full_tt(self, n_qubits):
        if n_qubits > 0:
            self.root = QubitNum(0)
            t = []
            level = 0
            qubits = list(range(0, n_qubits))
            while 3 ** level + (3 ** level - 1) // 2 < n_qubits:
                t = t + [(qubits[(3 ** level - 1) // 2:3 ** level + (3 ** level - 1) // 2])]
                level += 1

            t = t + [(qubits[(3 ** level - 1) // 2:n_qubits])]
            self.nodes = t
            self.update_branchnum()
            self[QubitNum(n_qubits - 1)][2] = LostNum()
        else:
            self._nodes = {}

    def build_st_tree(self, n_qubits):
        if n_qubits > 0:
            self.root = QubitNum(0)
            self._nodes[self.root] = NodeContacts(parent=ROOT,
                                                  childs=[BranchNum(1), BranchNum(2), 1])
            for i in range(1, n_qubits):
                self.nodes[QubitNum(i)] = NodeContacts(parent=i-1,
                                                       childs=[BranchNum(2*i+1), BranchNum(2*i+2), i+1])
            self._nodes[QubitNum(n_qubits - 1)][2] = LostNum()
            self.update_branchnum()
            self[QubitNum(n_qubits - 1)][2] = LostNum()
            self.check_branch_numeration()
        else:
            self._nodes = {}

    def check_branch_numeration(self):
        bnum = []
        for node in self:
            for i, child in enumerate(self[node]):
                if isinstance(child, BranchNum):
                    bnum.append(child.num)
        if set(bnum) != set(list(range(1, len(bnum) + 1))):
            raise TreeStructureError("Numeration of the branches should be from 1 to " + str(len(bnum)) + " without repetitions")


    def add_node(self, node, its_parent, pos = 0):

        if len(self.nodes) == 0:

            node = QubitNum(node)
            self.root = node
            self[node] = NodeContacts(ROOT)
        else:
            node = QubitNum(node)
            its_parent = QubitNum(its_parent)
            if node not in self.nodes and its_parent in self.nodes:
                if isinstance(self[its_parent][pos], (BranchNum, LostNum)):
                    self[node] = NodeContacts(its_parent)
                    self[its_parent][pos] = node
                else:
                    raise KeyError("Parent Node " + str(its_parent) +
                                   " already has child on " + str(pos))
            else:
                raise KeyError("Self.nodes already has node " + str(node) +
                               " or hasn't its_parent " + str(its_parent))

    def delete_node(self, nodes: list[int | QubitNum] | int | QubitNum = None):
        """
        Remove nodes and its childs from parent child
        """
        if isinstance(nodes, (QubitNum, int)):
            self._delete_node(nodes)
        else:
            for node in nodes:
                self._delete_node(node)

    def _delete_node(self, node: int | QubitNum = 0):
        """
        Remove node and its childs from nodes
        """
        if not isinstance(node, (QubitNum, int)):
            raise TypeError("Node should be integer or QubitNum")

        def down(
                node_: QubitNum
        ):
            for child_ in self[node_].childs:
                if not child_.is_last:
                    down(child_)
            self.nodes.pop(node_)

        node = QubitNum(node)
        if node in self.nodes:
            if self[node].parent == ROOT:
                self._nodes = {}
            else:
                for i, child in enumerate(self[self[node].parent]):
                    if node == child:
                        self.parent_contacts(node)[i] = LostNum()
            down(node)

    def parent_contacts(self, node):
        return self[self[node].parent]

    def update_branchnum(self, enum_list: list = None, renumerate = True):
        pos = 0

        def down(node):
            nonlocal pos, enum_list
            for i, child in enumerate(self[node]):
                if isinstance(child, BranchNum):
                    self[node][i] = BranchNum(enum_list[pos])
                    pos += 1
                elif not child.is_last:
                    down(child)

        if enum_list is None:
            pos = 1
            if renumerate:
                for node in self:
                    for i, child in enumerate(self[node]):
                        if isinstance(child, BranchNum):
                            self[node][i] = LostNum()
            for node in self:
                for i, child in enumerate(self[node]):
                    if isinstance(child, LostNum):
                        self[node][i] = BranchNum(pos)
                        pos += 1
        else:
            pos = 0
            try:
                down(self.root)
            except IndexError:
                raise TreeStructureError("Numeration should be from 1 to " + str(len(self.branches())))
            self.check_branch_numeration()

    def __getitem__(self, key):
        return self.nodes[QubitNum(key)]

    def __setitem__(self, key, val):
        self.nodes[QubitNum(key)] = val
        return val

    def __iter__(self):
        return self.nodes.__iter__()

    def branches(self) -> list[tuple[int, list[tuple[int, str]]]]:
        """
        Convert nodes to list of branches. Used by pauli tables
        """
        s = []

        def down(node, k, branches):
            for index, child in enumerate(self[node]):
                if not child.is_last:
                    down(child, k + [tuple([node.num, gate_number_to_name[index]])], branches)
                elif isinstance(child, BranchNum):
                    branches.append((child.num, k + [(node.num, gate_number_to_name[index])]))
        down(self.root, k=[], branches=s)
        return s

    #TODO
    def draw(self):
        pass

    def __str__(self, until_enum=False):
        k = []
        s = []
        enum_list = []
        L = self.height()

        def down(parent, k):
            nonlocal s, enum_list
            for index, child in enumerate(self[parent].childs):
                if not child.is_last:
                    down(child, k + [str(parent) + gate_number_to_name[index]])
                elif not isinstance(child, LostNum):
                    if not until_enum:
                        s.append(k + [str(parent) + gate_number_to_name[index]])
                        enum_list.append(self.nodes[parent][index])
                    else:
                        s.append(k)

        down(self.root, k)
        pr = ''
        k = 4
        if not until_enum:
            for l in range(0, L):
                for branch in s:
                    if len(branch) > l:
                        pr += branch[l] + ' ' * (k - len(branch[l]))
                    else:
                        pr += " " * k
                pr = pr + '\n'
            for num in enum_list:
                pr += str(num) + ' ' * (k - len(str(num)))
        else:
            for l in range(0, L - 1):
                for branch in s:
                    if len(branch) > l:
                        pr += branch[l] + ' ' * (k - len(branch[l]))
                        pass
                    else:
                        pr += " " * k
                pr = pr + '\n'

        return pr + '\n ---------------------------------------------------------------'

