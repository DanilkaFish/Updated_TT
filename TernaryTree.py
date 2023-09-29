from BaseTree import BaseTernaryTree,TreeStructureError
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

        def check_xy():
            """
            find unoccupied branch for transposition with Z...Z branch
            """
            for node in self.nodes:
                for index, child in enumerate(self[node]):
                    if not child:
                        return [node, index]
            raise TreeStructureError("Here is now lost branch")

        def descent_before_not_z(_node, _edge):
            _edge = (_edge + 1) % 2
            if not self[_node][_edge].is_last:
                _node = self[_node][_edge]
                _edge = 2
                while not self[_node][_edge].is_last:
                    _node = self[_node][_edge]
            return _node, _edge

        def rise_before_not_z(_node, _edge):
            while _edge == 2:
                # Подъем до не Z edge
                n = _node
                _node = self.parent(_node)
                _edge = self[_node].childs.index(n)
            return _node, _edge

        s = []
        # find Z...Z branch to eliminate it
        child = self[self.root][2]
        while not child.is_last:
            parent = child
            child = self[child][2]
        if child:
            first, index = check_xy()
            s.append(self.branch_transposition(first, index, parent, 2))

        nnodes = self.num_nodes

        for i in range(nnodes):
            branches = self.branches()
            try:
                num = {branches[j][0]: j for j in range(2*nnodes)}
            # print(num)
            except IndexError:
                raise IndexError("You should use 2n branches for numeration, where n = ",
                      self.num_nodes, "but you use",  len(branches))
            except:
                raise Exception("What's the hell??")

            node1, edge1 = branches[num[2*i + 1]][1][-1]
            node2, edge2 = branches[num[2*i + 2]][1][-1]
            edge1 = gate_index[edge1]
            edge2 = gate_index[edge2]

            node1, edge1 = rise_before_not_z(node1, edge1)
            node1, edge1 = descent_before_not_z(node1, edge1)

            r = self.branch_transposition(node1, edge1, node2, edge2)
            if len(r) > 0:
                s.append(r)
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
