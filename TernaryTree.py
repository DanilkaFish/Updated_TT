from BaseTree import BaseTernaryTree
from Nodes import QubitNum, BranchNum, LostNum

gate_name = {0: 'X', 1: 'Y', 2: 'Z'}
gate_index = {'X': 0, 'Y': 1, 'Z': 2}
dict_prod = {"II": "I", "XX": "I", "YY": "I", "ZZ": "I",
             "XY": 'Z', "YX": 'Z', "XZ": 'Y', "ZX": 'Y', "YZ": 'X', "ZY": 'X',
             "IX": "X", "XI": "X", "YI": "Y", "IY": "Y", "IZ": "Z", "ZI": "Z"}


class TernaryTree(BaseTernaryTree):
    """
    Ternary tree object for initial state preparation
    """
        
    def branch_transposition(
            self,
            node1: QubitNum | int,
            edge1: int,
            node2: QubitNum | int,
            edge2: int
    ):
        """
        Transpose subtrees or branches  wtih preserving of nodes numeration. Return pauli operator implementing this transformation.
        node1 = int -- node number in tree
        edge1 = 0|1|2 -- edge which connected node1 with its parent
        node2 = int -- node number in tree
        edge2 = 0|1|2 -- edge which connected node2 with its parent
        """
        ordered_qubits = {}
        qubits = sorted(self.nodes.keys())
        for i in range(self.num_nodes):
            ordered_qubits[qubits[i]] = i
        s = ["I"] * self.num_nodes

        node1 = QubitNum(node1)
        node2 = QubitNum(node2)
        parent_closest_node, _ = self._closest_parent(node1, node2)

        while parent_closest_node != node1:
            s[ordered_qubits[node1]] = gate_name[edge1]
            for index, child in enumerate(self[self.parent(node1)]):
                if child == node1:
                    edge1 = index
            node1 = self[node1].parent

        while parent_closest_node != node2:
            s[ordered_qubits[node2]] = gate_name[edge2]
            for index, child in enumerate(self[self.parent(node2)]):
                if child == node2:
                    edge2 = index
            node2 = self.parent(node2)
        s[ordered_qubits[parent_closest_node]] = dict_prod[gate_name[edge1] + gate_name[edge2]]
        
        # Tree transposition
        aux_node = self[node1][edge1]
        self[node1][edge1] = self[node2][edge2]
        if not self[node2][edge2].is_last:
            self[self[node2][edge2]].parent = node1

        self[node2][edge2] = aux_node
        if aux_node:
            self[aux_node].parent = node2
        return ''.join(s)

    def _closest_parent(self, node1: QubitNum | int, node2: QubitNum | int):
        """
        Return closest parent to the first and second node
        """
        def check_parent_child_relation(parent_node: QubitNum,
                       child_node: QubitNum
                       ):
            """
            Check whether init_node is parent of search_node or not
            """

            flag = True
            dist = 0

            def down(node_p: QubitNum | int,
                     node_c: QubitNum | int,
                     _dist: int ):
                nonlocal flag, dist
                for child in self[node_p].childs:
                    if not child.is_last:
                        if child == node_c:
                            dist = _dist + 1
                            flag = False

                            down(child, node_c, _dist + 1)

            if parent_node == child_node:
                return False, 0
            down(parent_node,child_node,0)
            return flag, dist

        flag, dist0 = check_parent_child_relation(node1, node2)
        i = 0
        while flag:
            node1 = self[node1].parent
            flag, dist0 = check_parent_child_relation(node1, node2)
            i += 1
        return node1, i + dist0

    def to0vac(self):
        """
        Use algorithm described in the graduate work
        """
        help_dict = {0:1, 1:0}

        def _find_second_branch(node,edge):
            """
            This function search branch to tranpose with given 
            """
            parent2, flag = -1, False

            def down(parent):
                nonlocal parent2, flag
                for index, child in enumerate(self.nodes[parent].childs):
                    if flag:
                        break
                    if parent == node and index == edge:
                        i=0
                    elif isinstance(child,EnumInfo) and index == 2:                        
                        parent2 = parent
                        flag = True
                    elif child:
                        down(child)
            down(0)
            return parent2

        def _check_xy():
            """
            find unoccupied branch for transposition with Z...Z branch  
            """
            parentf = 0 
            indexf = 0
            def down(parent):
                nonlocal parentf,indexf 
                for index, child in enumerate(self.nodes[parent].childs):
                    if isinstance(child, bool):
                        parentf, indexf = parent, index 
                    if child:
                        down(child)
            down(0)
            return [parentf,indexf]
        
        def check_pair(pair):
            if min(pair[0],pair[1]) % 2 == 1 and abs(pair[0] - pair[1]) == 1:
                return True
            else:
                return False

        def find_trans(pair_list,index):
            pair = pair_list[index]
            indexes = []
            if check_pair(pair):
                return False
            for i, ipair in enumerate(pair_list[index + 1:]):
                if check_pair([ipair[0],pair[0]]):
                    return [2*index + 1, 2*(i + index + 1) ]
                if check_pair([ipair[0],pair[1]]):
                    return [2*index, 2*(i + index + 1) ]
                if check_pair([ipair[1],pair[0]]):
                    return [2*index + 1, 2*(i + index + 1) + 1]
                if check_pair([ipair[1],pair[1]]):
                    return [2*index, 2*(i + index + 1) + 1]
       
                
        def get_pair(self,node1,edge1, node2, edge2, branches):
            l = [branch[-1] for branch in branches]
            index = self.enum_list[l.index([node1,gate_name[edge1]])]
            if index % 2 ==0:
                _index = index + 1
            else:
                _index = index - 1
            _node, _edge = branches[self.enum_list.index(_index)][-1]
            _, dist = self._closest_parent(node2,_node) 
            return _node, gate_index[_edge], dist
        
        
        
        s = []
        parent = 0 
        
        # find Z...Z branch to eliminate it
        child = self.nodes[parent][2]
        while child:
            parent = child
            child = self.nodes[parent][2]
        if not isinstance(child,bool):
            first, index = _check_xy()
            s.append(self.branch_transposition(first, index, parent,2))
        
#         Лучше через nodes сделать
#         first leaf
        i = 0
        for i in range(self.nmodes):
            branches = self.branches()
            node, edge = branches[self.enum_list.index(2*i)][-1]
            edge = gate_index[edge]
            if edge == 0:
                if isinstance(self.nodes[node][1],bool):
                    s.append(self.branch_transposition(node,1,node,2))
                if isinstance(self.nodes[node][1],EnumInfo):
                    _node = node
                    _edge = 1
                
                else:
                    _node = self.nodes[node][1]
                    _edge = 2
                    while not isinstance(self.nodes[_node][2],EnumInfo):
                        if isinstance(self.nodes[_node][2],bool):
                            s.append(self.branch_transposition(_node,2,_node,1))
                        else:
                            _node = self.nodes[_node][2]
                branches = self.branches()
                
                node1,edge1,dist1 = get_pair(self,node, edge, _node, _edge, branches)
#                 node2,edge2,dist2 = get_pair(self,_node, _edge, node, edge, branches) 
                if node1 != _node or edge1 != _edge:
                    
#                     if dist1 < dist2:
                    s.append(self.branch_transposition(_node,_edge,node1,edge1))
#                     else:
#                         s.append(self.branch_transposition(node,edge,node2,edge2))
            if edge == 1:
                if isinstance(self.nodes[node][0],bool):
                    s.append(self.branch_transposition(_node,0,node,2))
                if isinstance(self.nodes[node][0],EnumInfo):
                    _node = node
                    _edge = 0
                
                else:
                    _node = self.nodes[node][0]
                    _edge = 2
                    while not isinstance(self.nodes[_node][2],EnumInfo):
                        if isinstance(self.nodes[node][2],bool):
                            s.append(self.branch_transposition(_node,2,_node,0))
                        else:
                            _node = self.nodes[_node][2]
                branches = self.branches()
                
                node1,edge1,dist1 = get_pair(self,node, edge, _node, _edge, branches)
#                 node2,edge2,dist2 = get_pair(self,_node, _edge, node, edge, branches) 
                if node1 != _node or edge1 != _edge:
                    
#                     if dist1 < dist2:
                    s.append(self.branch_transposition(_node,_edge,node1,edge1))
#                     else:
#                         s.append(self.branch_transposition(node,edge,node2,edge2))
            if edge == 2:
                n = node
                _node = self.nodes[node].parent
                _edge = self.nodes[_node].childs.index(node)
                while _edge == 2:
                    # Подъем до не Z edge
                    n = _node
                    _node = self.nodes[_node].parent
                    _edge = sefl.nodes[_node].childs.index(n)
                n = _node
                # cпуск по Z
                if _edge == 0:
                    if isinstance(self.nodes[n][1],bool):
                        s.append(self.branch_transposition(n,1,n,2))
                        
                    if isinstance(self.nodes[n][1],EnumInfo):
                        _node = n
                        _edge = 1

                    else:
                        _node = self.nodes[n][1]
                        _edge = 2
                        while not isinstance(self.nodes[_node][2],EnumInfo):
                            if isinstance(self.nodes[_node][2],bool):
                                s.append(self.branch_transposition(_node,2,_node,1))
                            else:
                                _node = self.nodes[_node][2]
                    
                    branches = self.branches()
                    node1,edge1,dist1 = get_pair(self,node, edge, _node, _edge, branches)
                    node2,edge2,dist2 = get_pair(self,_node, _edge, node, edge, branches) 
                    if node1 != _node or edge1 != _edge:

                        if dist1 < dist2:
                            s.append(self.branch_transposition(_node,_edge,node1,edge1))
                        else:
                            s.append(self.branch_transposition(node,edge,node2,edge2))
                if edge == 1:
                    if isinstance(self.nodes[n][0],bool):
                        s.append(self.branch_transposition(n,0,n,2))
                    if isinstance(self.nodes[n][0],EnumInfo):
                        _node = node
                        _edge = 0

                    else:
                        _node = self.nodes[n][0]
                        _edge = 2
                        while not isinstance(self.nodes[_node][2],EnumInfo):
                            if isinstance(self.nodes[_node][2],bool):
                                s.append(self.branch_transposition(_node,2,_node,0))
                            else:
                                _node = self.nodes[_node][2]
                    branches = self.branches()
                    node1,edge1,dist1 = get_pair(self,node, edge, _node, _edge, branches)
                    node2,edge2,dist2 = get_pair(self,_node, _edge, node, edge, branches) 
                    if node1 != _node or edge1 != _edge:

                        if dist1 < dist2:
                            s.append(self.branch_transposition(_node,_edge,node1,edge1))
                        else:
                            s.append(self.branch_transposition(node,edge,node2,edge2))
        return s
    
    def tojw(self):
        """
        Use less effective algorithms with transformation to JW tree
        """
        s = []
        def _find_first_branch(node):
            if self.nodes[node][0]:
                return [node,0]
            if self.nodes[node][1]:
                return [node,1]
            if self.nodes[node][2]:
                return _find_first_branch(self.nodes[node][2])
            return [False, False]
    # ищет первую не занятую z вершину   
        def _find_second_branch(node,edge):
            parents,flag = -1, False
            def down(parent):
                nonlocal parents, flag
                for index, child in enumerate(self.nodes[parent].childs):
                    if flag:
                        break
                    if parent == node and index == edge:
                        i=0
                    elif isinstance(child,EnumInfo) and index == 2:                        
                        parents = parent
                        flag = True
                    elif child:
                        down(child)
            down(0)
            return parents
        def _check_xy():
            parentf = 0 
            indexf = 0
            def down(parent):
                nonlocal parentf,indexf 
                for index, child in enumerate(self.nodes[parent].childs):
                    if isinstance(child, bool):
                        parentf, indexf = parent, index 
                    if child:
                        down(child)
            down(0)
            return [parentf,indexf]
        def check_pair(pair):
            if min(pair[0],pair[1]) % 2 == 1 and abs(pair[0] - pair[1]) == 1:
                return True
            else:
                return False
        def find_trans(pair_list,index):
            pair = pair_list[index]
            indexes = []
            if check_pair(pair):
                return False
            for i, ipair in enumerate(pair_list[index + 1:]):
                if check_pair([ipair[0],pair[0]]):
                    return [2*index + 1, 2*(i + index + 1) ]
                if check_pair([ipair[0],pair[1]]):
                    return [2*index, 2*(i + index + 1) ]
                if check_pair([ipair[1],pair[0]]):
                    return [2*index + 1, 2*(i + index + 1) + 1]
                if check_pair([ipair[1],pair[1]]):
                    return [2*index, 2*(i + index + 1) + 1]
                
        first, index = _check_xy()
#       Здесь сводится все к дереву JW
        s.append(self.branch_transposition(first, index, first, 2))
        s.append(self.branch_transposition(first, 0, _find_second_branch(first,0), 2))
        s.append(self.branch_transposition(first, 1, _find_second_branch(first,1), 2))
        while True:
            fl = _find_first_branch(0)
            if fl[0] or not isinstance(fl[0], bool):
                s.append(self.branch_transposition(fl[0], fl[1], _find_second_branch(fl[0],fl[1]), 2))
            else:
                break
#       Далее просто переставляются оставшиеся ветви
        branches, num_list = self.branches(True)
        num_list = [i.num for i in num_list]
        pair_list = [[num_list[2*i], num_list[2*i+1]] for i in range(len(num_list)//2)]
        L = len(pair_list)
        for index in range(L):
            k = find_trans(pair_list, index)
            if k:
                q1 = branches[k[0]][-1]
                q2 = branches[k[1]][-1]
                s.append(self.branch_transposition(q1[0], gate_index[q1[1]], q2[0], gate_index[q2[1]]))
            branches, num_list = self.branches(True)
            num_list = [i.num for i in num_list]
            pair_list = [[num_list[2*i], num_list[2*i+1]] for i in range(len(num_list)//2)]
        return s


def prod_pauli_strings(pauli1: list[tuple[int, str]],
                       pauli2: list[tuple[int, str]]):
    pauli1 = {gate[0]: gate[1] for gate in pauli1}
    pauli2 = {gate[0]: gate[1] for gate in pauli2}
    pauli = {}
    for qubit in set(list(pauli1.keys()) + list(pauli2.keys())):
        pauli[qubit] = dict_prod[pauli1.get(qubit, "I") + pauli2.get(qubit, "I")]

    return list(pauli.items())


def pauli_weight(tt: TernaryTree, indexes: list[int]):
    branches = tt.branches()
    branches = [branch[1] for branch in branches if branch[0] in indexes]
    prod = branches[0]
    for pauli in branches[1:]:
        prod = prod_pauli_strings(prod, pauli)
    return len([p for p in prod if p[1] != "I"])
