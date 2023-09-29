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
        node1 = QubitNum -- node number in tree
        edge1 = 0|1|2 -- edge which connected node1 with its parent
        node2 = QubitNum -- node number in tree
        edge2 = 0|1|2 -- edge which connected node2 with its parent
        """
        ordered_qubits = {}
        qubits = sorted(self.nodes.keys())
        for i in range(self.num_nodes):
            ordered_qubits[qubits[i]] = i
        s = ["I"] * self.num_nodes
        _edge1 = edge1
        _edge2 = edge2
        _node1 = QubitNum(node1)
        _node2 = QubitNum(node2)
        node1 = QubitNum(node1)
        node2 = QubitNum(node2)
        if node1 == node2 and edge1 == edge2:
            return ""

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
        aux_node = self[_node1][_edge1]
        self[_node1][_edge1] = self[_node2][_edge2]
        if not self[_node2][_edge2].is_last:
            self[self[_node2][_edge2]].parent = _node1
        self[_node2][_edge2] = aux_node
        if not aux_node.is_last:
            self[aux_node].parent = _node2
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
            parent = child_node
            dist = 0
            while not parent.num == self.root.num:
                if parent == parent_node:
                    return True, dist
                child_node = parent
                parent = self.parent(child_node)
                dist += 1
            if parent == parent_node:
                return True, dist
            return False, 0

        flag, dist0 = check_parent_child_relation(node1, node2)
        i = 0
        while not flag:
            node1 = self[node1].parent
            flag, dist0 = check_parent_child_relation(node1, node2)
            i += 1
        return node1, i + dist0

    def to0vac(self):
        """
        Use algorithm described in the graduate work
        """
        help_dict = {0:1, 1:0}

        # def _find_second_branch(node, edge):
        #     """
        #     This function search branch to tranpose with given
        #     """
        #     parent2, flag = -1, False
        #
        #     def down(_node):
        #         nonlocal parent2, flag
        #         for index, child in enumerate(self[_node].childs):
        #             if flag:
        #                 break
        #             if _node == node and index == edge:
        #                 i=0
        #             elif isinstance(child,BranchNum) and index == 2:
        #                 parent2 = _node
        #                 flag = True
        #             elif not child:
        #                 down(child)
        #     down(0)
        #     return parent2

        def _check_xy():
            """
            find unoccupied branch for transposition with Z...Z branch  
            """
            for node in self.nodes:
                for index, child in enumerate(self[node]):
                    if not child:
                        return [node, index]
        
        # def check_pair(pair):
        #     if min(pair[0], pair[1]) % 2 == 1 and abs(pair[0] - pair[1]) == 1:
        #         return True
        #     else:
        #         return False
        #
        # def find_trans(pair_list,index):
        #     pair = pair_list[index]
        #     indexes = []
        #     if check_pair(pair):
        #         return False
        #     for i, ipair in enumerate(pair_list[index + 1:]):
        #         if check_pair([ipair[0],pair[0]]):
        #             return [2*index + 1, 2*(i + index + 1) ]
        #         if check_pair([ipair[0],pair[1]]):
        #             return [2*index, 2*(i + index + 1) ]
        #         if check_pair([ipair[1],pair[0]]):
        #             return [2*index + 1, 2*(i + index + 1) + 1]
        #         if check_pair([ipair[1],pair[1]]):
        #             return [2*index, 2*(i + index + 1) + 1]
        #
        # def get_pair(self, node1, edge1, node2, edge2, branches):
        #     l = [branch[-1] for branch in branches]
        #     index = self.enum_list[l.index([node1,gate_name[edge1]])]
        #     if index % 2 ==0:
        #         _index = index + 1
        #     else:
        #         _index = index - 1
        #     _node, _edge = branches[self.enum_list.index(_index)][-1]
        #     _, dist = self._closest_parent(node2,_node)
        #     return _node, gate_index[_edge], dist
        
        
        
        s = []

        # find Z...Z branch to eliminate it
        child = self[self.root][2]
        while not child.is_last:
            parent = child
            child = self[child][2]
        if not child:
            first, index = _check_xy()
            s.append(self.branch_transposition(first, index, parent, 2))

        # Лучше через nodes сделать
        # first leaf
        nnodes = self.num_nodes
        for i in range(nnodes):
            branches = self.branches()
            num = {branches[j][0]: j for j in range(2*nnodes)}
            node1, edge1 = branches[num[2*i + 1]][1][-1]
            node2, edge2 = branches[num[2*i + 2]][1][-1]
            edge1 = gate_index[edge1]
            edge2 = gate_index[edge2]
            if edge1 == 0:
                if self[node1][1].is_last:
                    r = self.branch_transposition(node1, 1, node2, edge2)
                    if len(r) > 0:
                        s.append(r)
                else:
                    _node = self[node1][1]

                    while not self[_node][2].is_last:
                        _node = self[_node][2]
                    r = self.branch_transposition(_node, 2, node2, edge2)
                    if len(r) > 0:
                        s.append(r)

            if edge1 == 1:
                if self[node1][0].is_last:
                    r = self.branch_transposition(node1, 0, node2, edge2)
                    if len(r) > 0:
                        s.append(r)
                else:
                    _node = self[node1][0]

                    while not self[_node][2].is_last:
                        _node = self[_node][2]
                    r = self.branch_transposition(_node, 2, node2, edge2)
                    if len(r) > 0:
                        s.append(r)

            if edge1 == 2:
                n = node1
                node1 = self.parent(node1)
                edge1 = self[node1].childs.index(n)
                while edge1 == 2:
                    # Подъем до не Z edge
                    n = node1
                    node1 = self.parent(node1)
                    edge1 = self[node1].childs.index(n)

                # cпуск по Z
                if edge1 == 0:
                    if self[node1][1].is_last:
                        r = self.branch_transposition(node1, 1, node2, edge2)
                        if len(r) > 0:
                            s.append(r)
                    else:
                        _node = self[node1][1]

                        while not self[_node][2].is_last:
                            _node = self[_node][2]
                        r = self.branch_transposition(_node, 2, node2, edge2)
                        if len(r) > 0:
                            s.append(r)

                if edge1 == 1:
                    if self[node1][0].is_last:
                        r = self.branch_transposition(node1, 0, node2, edge2)
                        if len(r) > 0:
                            s.append(r)
                    else:
                        _node = self[node1][0]

                        while not self[_node][2].is_last:
                            _node = self[_node][2]
                        r = self.branch_transposition(_node, 2, node2, edge2)
                        if len(r) > 0:
                            s.append(r)
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
                    elif isinstance(child,BranchNum) and index == 2:                        
                        parents = parent
                        flag = True
                    elif not child:
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
                    if not child:
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
    d = {}
    for i in indexes:
        d[i] = d.get(i,0) + 1

    indexes = [key for key in d if d[key] % 2 == 1]
    if len(indexes) > 0:
        branches = [branch[1] for branch in branches if branch[0] in indexes]
        prod = branches[0]
        for pauli in branches[1:]:
            prod = prod_pauli_strings(prod, pauli)
        # print(indexes,":  ", len([p for p in prod if p[1] != "I"]))
        return len([p for p in prod if p[1] != "I"])
    return 0
