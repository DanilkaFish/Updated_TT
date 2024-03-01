from typing import Dict, List

from qiskit.circuit import Parameter, QuantumCircuit, ParameterExpression, QuantumRegister
from numpy import pi
from qiskit.circuit.library import HamiltonianGate
from qiskit.quantum_info import Operator


class MyCirq(QuantumCircuit):

    def __init__(self, n):
        self.paulis = []
        self.swaps = []
        # qr = QuantumRegister(n)
        super().__init__(n)

    def pauli(self, pl: Dict[int, str], qubs: List[int], par: float | Parameter = - pi/4, sign=1):
        if not isinstance(pl, dict):
            pl = {gate[0]: gate[1] for gate in pl}
        pauli = {}
        for el in pl:
            if pl[el] != "I":
                pauli[el] = pl[el]

        label = ''.join(pauli.values())
        op = Operator.from_label(label)
        qubits = list([i for i in pauli.keys()])
        gate = HamiltonianGate(op, par*sign)
        if isinstance(par, float):
            # gate.name = 'pi/4'
            # "Меня зовут %s. Мне %d лет." % (name, age)
            gate.name = "exp{i %s pi/4}" %label
            pass
        if isinstance(par, Parameter):
            # gate.name = str(par)[6:10]
            gate.name = '0'
            pass
        if isinstance(par, ParameterExpression):
            gate.name = str(par)
            gate.name = '0'
            pass
        self.append(gate, qubits)
        self.paulis.append((pauli, par))

