from __future__ import annotations


from Nodes import BranchNum, QubitNum, LostNum, NodeContacts, R


gate_name_to_number = {'X': 0, 'Y': 1, 'Z': 2}
gate_number_to_name = {0: 'X', 1: 'Y', 2: 'Z'}
dict_prod = {"II": "I", "XX": "I", "YY": "I", "ZZ": "I",
             "XY": 'Z', "YX": 'Z', "XZ": 'Y', "ZX": 'Y', "YZ": 'X', "ZY": 'X',
             "IX": "X", "XI": "X", "YI": "Y", "IY": "Y", "IZ": "Z", "ZI": "Z"}


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

    def is_root(self, node: QubitNum | int):
        if self[node].parent == R:
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

    # TODO
    @nodes.setter
    def nodes(self, nodes):

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

    def check_branch_numeration(self):
        bnum = []
        for node in self:
            for i, child in enumerate(self[node]):
                if isinstance(child, BranchNum):
                    bnum.append(child.num)
        if set(bnum) != set(list(range(1, len(bnum) + 1))):
            raise ValueError("Numeration of the branches should be from 1 to " + str(len(bnum)) + " without repetitions")

    def check_qubits_structure(self):
        """
            Small check of nodes object.
        """

        # check = {}
        # edges = {}

        def _check_qubits_structure(child, **kwargs):
            nodes = kwargs["nodes"]
            parent = kwargs["parent"]

            if isinstance(child, QubitNum):
                if not self.is_root(child):
                    return nodes[child].parent == parent

        s = self.recursive_traversal(_check_qubits_structure)
        for obj in s:
            if isinstance(obj, bool):
                if not obj:
                    raise ValueError("Some of the child's parent and parent's child differ")
        s = [key for key in self.nodes]
        if len(set(s)) != len(s):
            raise ValueError("Numeration of the qubits should differ")

    def add_node(self, node, its_parent, pos = 0):

        if len(self.nodes) == 0:

            node = QubitNum(node)
            self.root = node
            self[node] = NodeContacts(R)
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

    def delete_nodes(self, nodes: list[int|QubitNum] = None):
        """
        Remove nodes and its childs from parent child
        """
        for node in nodes:
            self.delete_node(node)

    def delete_node(self, node=0):
        """
        Remove node and its childs from nodes
        """
        def down(
                node_: QubitNum
        ):
            for child_ in self[node_].childs:
                if not child_.is_last:
                    down(child_)
            self.nodes.pop(node_)

        node = QubitNum(node)
        if node in self.nodes:
            if self[node].parent == R:
                self.nodes = {}
            else:
                for i, child in enumerate(self[self[node].parent]):
                    if node == child:
                        self.parent_contacts(node)[i] = LostNum()
            down(node)
        # else:
        #     warnings.warn("There is already no node " + str(node) + " in self.nodes")




    def parent_contacts(self, node):
        return self[self[node].parent]

    def build_st_tree(self, n_qubits):
        self._nodes[self.root] = NodeContacts(parent=R,
                                              childs=[BranchNum(1), BranchNum(2), 1])
        for i in range(1, n_qubits):
            self.nodes[QubitNum(i)] = NodeContacts(parent=i-1,
                                                   childs=[BranchNum(2*i+1), BranchNum(2*i+2), i+1])
        self._nodes[QubitNum(n_qubits - 1)][2] = LostNum()
        self.check_branch_numeration()
        self.check_qubits_structure()


    def update_branchnum(self, enum_list: list = None):
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
            down(self.root)
            self.check_branch_numeration()

    def __getitem__(self, key):
        return self.nodes[QubitNum(key)]
    # def __setitem__(self, key, val):
    #     childs = list(self.childs)
    #     childs[key] = self.obj_to_child(val)
    #     self.childs = childs
    #     return val

    def __setitem__(self, key, val):
        self.nodes[QubitNum(key)] = val
        return val

    def __iter__(self):
        return self.nodes.__iter__()

    # def __next__(self):
    #     return self.nodes.__next__()

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
        num_space = 0
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
        #             for num in enum_list:
        #                 pr += str(num) + ' '*(k-len(str(num)))
        return pr + '\n ---------------------------------------------------------------'

    def recursive_traversal(self, fun, init=None, **kwargs):
        if init is None:
            init = self.root
        s_ = []
        if 'nodes' in kwargs:
            nodes = kwargs["nodes"]
        else:
            nodes = self._nodes

        def down(
                s: list[int],
                parent: QubitNum
        ):
            for i, child in enumerate(nodes[parent]):
                s.append(fun(child, parent=parent, i=i, nodes=nodes, **kwargs))
                if not child.is_last:
                    down(s, child)
        try:
            down(s_, QubitNum(init))
        except KeyError:
            raise KeyError("nodes dictionary should content all the child nodes and root node")
        except Exception as e:
            raise e
        return s_
    # @property
    # def enum_list(self):
    #     """
    #     Branches numeration
    #     """
    #     return self._enum_list
    
    # @enum_list.setter
    # def enum_list(self,num_list):
    #     def check_enum_list():
    #         inf = min(self.enum_list)
    #         if len(self._enum_list) == self.nmodes*2 :
    #             self._enum_list = [num - inf for num in self.enum_list]
    #         else:
    #             raise ValueError("Wrong numeration list length. Should be " +
    #             str(2*self.nmodes) + " , but " + str(len(self._enum_list)) + " were given ")
    #         for i in range(2*self.nmodes):
    #             if i not in self.enum_list:
    #                 raise ValueError("Numeration has to contain all the number beetwen 0 and " +
    #                 str(2*self.nmodes) + ", but yout list is " + str(self.enum_list))
    #     if num_list:
    #         self._enum_list = num_list
    #     else:
    #         self._enum_list = st_enumeration(self.nmodes)
    #     if isinstance(self._enum_list[0], list) :
    #         l = []
    #         for pair in self.enum_list:
    #             l += pair
    #         self._enum_list = l
    #     check_enum_list()
    #     self.num_branches()
        
    # @property
    # def min_height(self):
    #     """
    #     Return minimal possible tree's height for given number of qubits
    #     """
    #     return int(trunc(log(2*self.n_qubits + 1)/log(3) - 0.001)) + 1
    #
    #
    
#     def _check_nodes(self,nodes):
#         "Small check of nodes object"
# #             nodes = self.nodes
#         check = {}
#         edges = {}
#         def append(parent):
#             if not isinstance(parent, int):
#                 raise ValueError("Parent should be int object, not " + str(type(parent)))
#             if parent in check:
#                 raise ValueError('Multiple initialization of parent node with number ' + str(parent))
#             check[parent] = None
#
#         for parent in nodes:
#             append(parent)
#             if isinstance(nodes[parent], NodeContacts):
#                 pass
#             else:
#                 raise ValueError("KeyValue in self.nodes should be NodeContacts object, not " + str(type( nodes[parent])))
#
        
#     def set_data(self, user_change = True):
#         """
#         Set n_qubits and nmodes from nodes info
#         """
#         self.n_qubits = len(self.nodes)
#         self.nmodes = 0
#         num_list = []
#         renum_flag = False
#         parent = 0
#         def down(parent):
#             nonlocal renum_flag
#             for child in self.nodes[parent].childs:
#                 if isinstance(child,EnumInfo):
#                     num_list.append(child.num)
#                     self.nmodes += 1
#                 elif child:
#                     down(child)
#                 elif child is None:
#                     renum_flag = True
#                     self.nmodes += 1
#         down(parent)
# #         if self.nmodes % 2 == 1:
# #             raise ValueError("nodes should contain only even number branches")
#         self.nmodes = self.nmodes // 2
#         if renum_flag or not user_change:
#
#             self.enum_list = st_enumeration(self.nmodes)
#         else:
#
#             self.enum_list = num_list
            
#     def check_height(self):
#         """
#         This method should be wherever the height could be changed, Set BaseTree._height
#         """
#         h = 0
#         h_max = 0
#         def down(parent,h):
#             h+= 1
#             nonlocal h_max
#             for index, child in enumerate(self.nodes[parent].childs):
#                 if child:
#                     down(child, h)
#                 elif h > h_max:
#                     h_max = h
#         down(0,h)
#         self._height = h_max
#
#
#     def set_num_qubit_transform(self):
#         """
#         Function for qubit consistent numeration due to its possible deleting or inserting.
#         """
#         self.num_qubit_transform = {}
#         for i, parent in enumerate(self.nodes):
#             self.num_qubit_transform[parent] = i
#
#     def renum_nodes(self):
#         """
#         Renumerate nodes in nodes after after changing the tree
#         """
#         self.set_num_qubit_transform()
#         for parent in list(self.nodes):
#             # Changing KeyValue
#             for index, child in enumerate(self.nodes[parent].childs):
#                 if child:
#                     self.nodes[parent][index] = self.num_qubit_transform[child]
#             if self.nodes[parent].parent >= 0:
#                 self.nodes[parent].parent = self.num_qubit_transform[self.nodes[parent].parent]
#             # Creating new renumerated key and deleting old key
#             self.nodes[self.num_qubit_transform[parent]] = self.nodes.pop(parent)
#         self.set_num_qubit_transform()
#
#     # Old method for tree inizialization
#     def build_alpha_beta_tree(self, height = 0):
#         r"""
#             Supports only even number electrons. Divide tree on 2 equals parts for alpha and beta electrons.
#         """
#         def _delete_node(self,node):
#             """
#             Remove nodes and its childs from parent child
#             """
#             def erase(parent):
#                 if parent:
#                     for child in nodes[parent].childs:
#                         if child:
#                             erase(child)
#                     nodes.pop(parent)
#
#             def erase_from_parent(child):
#                 parent = nodes[child].parent
#                 if parent in nodes:
#                     i = nodes[parent].childs.index(child)
#                     nodes[parent].childs[i] = None
#             if node in nodes:
#                 erase_from_parent(node)
#                 erase(node)
#
#         nodes = self._build_max_tree(height)
#         num_nodes = len(nodes) - self.n_qubits
#         L = self.min_height
#         l = L
#         center = (3**(L-1) - 1)//2 + 3**(L-1)//2
#         parent = 0
#         while (3**l - 1)//2 > num_nodes:
#             parent = nodes[parent][1]
#             l -= 1
#         left = center - 3**(l-1)//2 - 1
#         right = center + 3**(l-1)//2 + 1
#
#         if ((3**l - 1)//2 - num_nodes) %2 == 0:
#             _delete_node(self,parent)
#             num_nodes = num_nodes - (3**l - 1)//2
#         else:
#             for child in nodes[parent]:
#                 _delete_node(self,child)
#             num_nodes = num_nodes - (3**l - 1)//2 + 1
#
#         i = 0
#         while len(nodes) > self.n_qubits:
#             _delete_node(self,right + i)
#             _delete_node(self,left - i)
#             i += 1
#         parent = 0
#         while nodes[parent][1]:
#             parent = nodes[parent][1]
#         nodes[parent][1] = False
#         self.nodes = nodes
#         self.set_num_qubit_transform()
#         self.renum_nodes()
#
#     def _build_max_tree(self, height = 0):
#         """
#         Build tree (nodes) with all possible nodes and childs with height = self.min_height | height. Сonvenient for obtaining the necessary structures through nodes removal.
#         """
#         nodes = {}
#         L = max(height,self.min_height)
#         full_nodes = (3**L - 1) // 2
#         n = 0
#         first_parent = 0
#         first_child = 1
#         for l in range(L - 1):
#             for parent in range(first_parent, first_parent + 3**l):
#                 nodes[parent] = NodeContacts(first_parent + (parent - first_parent ) // 3 - 3**(l-1))
#                 for child in range(3):
#                     nodes[parent][child] = first_child + child
#                 first_child += 3
#             first_parent = first_parent + 3**l
#         for parent in range(first_parent, first_parent + 3**(L-1)):
#             nodes[parent] = NodeContacts(first_parent + (parent - first_parent ) // 3 - 3**(L-2))
#         return nodes
#
#
#     @classmethod
#     def build_max_tree(cls, height = 1, num_false = None,edge_false = None):
#         """
#         Build tree with all possible nodes and childs with height = self.min_height | height. Сonvenient for obtaining the necessary structures through nodes removal.
#         """
#         nodes = {}
#         L = height
#         full_nodes = (3**L - 1) // 2
#         if num_false is None:
#             num_false = full_nodes - 1
#         if edge_false is None:
#             edge_false = 2
#         if num_false not in range(full_nodes - 3**(L-1), full_nodes) or edge_false not in [0,1,2]:
#             raise ValueError("num_false and edge_false should to be in [" + str(full_nodes - 3**(L-1)) + ",...," + str(full_nodes -1) + "] and [0,1,2] respectively")
#         n = 0
#         first_parent = 0
#         first_child = 1
#         for l in range(L - 1):
#             for parent in range(first_parent, first_parent + 3**l):
#                 nodes[parent] = NodeContacts(first_parent + (parent - first_parent ) // 3 - 3**(l-1))
#                 for child in range(3):
#                     nodes[parent][child] = first_child + child
#                 first_child += 3
#             first_parent = first_parent + 3**l
#         for parent in range(first_parent, first_parent + 3**(L-1)):
#             nodes[parent] = NodeContacts(first_parent + (parent - first_parent ) // 3 - 3**(L-2))
#         nodes[num_false][edge_false] = False
#         tt = cls(nodes = nodes)
#         return tt
#
#
#     def num_branches(self,enum_list = None):
#         """
#         Numerate branches of the tree
#         """
#         if enum_list:
#             self.enum_list = enum_list
#
#         k = []
#         s = []
#         i = 0
#         def down(parent,k):
#             nonlocal s, i
#             for index, child in enumerate(self.nodes[parent].childs):
#                 if child:
#                     down(child, k +  [[parent, gate_name[index]]])
#                 elif (child != False) :
#                     if (i < len(self.enum_list)):
#                         self.nodes[parent][index] = EnumInfo(self.enum_list[i] + 1)
#                         i += 1
#                     else:
#                         self.nodes[parent][index] = False
#         down(0,k)
#         return self
#
#     def add_node(self, number = 1, parent = 0, gate = "Z"):
#         if number in self.nodes:
#             raise ValueError("node with number" + str(number) + "already exists")
#         if parent not in self.nodes:
#             raise ValueError("parent node with number" + str(number) + "does't exists")
#         if isinstance(self.nodes[parent].childs[gate_number[gate]], (bool,EnumInfo)):
#             childs = [None]*3
#             self.nodes[number] =  NodeContacts(parent, childs )
#             self.nodes[parent].childs[gate_number[gate]] = number
#
#         else:
#             raise ValueError("place " + gate + " is zanat")
#
#         self.renum_nodes()
#         self.set_data(user_change = False)
#         self.check_height()
#         self.num_branches()
#
#     def delete_node(self,node = 0 , nodes = None):
#         """
#         Remove nodes and its childs from parent child
#         """
#         flag = False
#         def erase(parent):
#             nonlocal flag
#             if parent:
#                 for child in self.nodes[parent].childs:
#                     if child:
#                         erase(child)
#                     elif child == False:
#                         flag = True
#                 self.nodes.pop(parent)
#
#         def parent_ch_num(child):
#             parent = self.nodes[child].parent
#             if parent in self.nodes:
#                 i = self.nodes[parent].childs.index(child)
#                 self.nodes[parent].childs[i] = None
#                 return parent, i
#         if nodes is None:
#             nodes = [node]
#         for node in nodes:
#             if node in self.nodes:
#                 parent, i = parent_ch_num(node)
#                 erase(node)
#                 if flag:
#                     self.nodes[parent].childs[i] = False
#             else:
#                 raise ValueError("There is no node " + str(node) + " in tree")
#         self.renum_nodes()
#         self.set_data(user_change = False)
#         self.check_height()
#         self.num_branches()
# #         print(self)
    
    

    
   
    