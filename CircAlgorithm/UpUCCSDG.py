from __future__ import annotations

import itertools

from . import TernaryTree, prod_pauli_strings
from .AbstractUCC import AbstractUCC
from .Mycirq import MyCirq


class UpUCCSDG(AbstractUCC):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tt = TernaryTree(self.num_spin_orbitals)

    def get_alpha_excitations(self) -> list[tuple[int, int]]:
        """
        Method to get spin preserving single excitation operators via ladder
        operators for alpha orbitals
        """
        alpha_occ = range(self.num_alpha)
        alpha_unocc = range(self.num_alpha, self.num_spatial_orbitals)
        return list(itertools.product(alpha_unocc, alpha_occ))

    def get_beta_excitations(self) -> list[tuple[int, int]]:
        """
        Method to get spin preserving single excitation operators via ladder
        operators for beta orbitals
        """
        beta_index_offset = self.num_spatial_orbitals
        beta_occ = range(beta_index_offset, beta_index_offset + self.num_beta)
        beta_unocc = range(beta_index_offset + self.num_beta, self.num_spin_orbitals)
        return list(itertools.product(beta_unocc, beta_occ))

    def get_double_excitations(self):
        """
        Method to get spin preserving pair double excitation operators via ladder
        operators
        """
        alpha = range(self.num_spatial_orbitals)
        beta_index_offset = self.num_spatial_orbitals
        ab = list(itertools.product(alpha, alpha))
        ab = [(el[0] + beta_index_offset, el[0], el[1] + beta_index_offset, el[1]) for el in ab]
        return ab

    def get_parametrized_circuit(self):
        n = self.num_spin_orbitals
        cirq = MyCirq(n)

        qubs = [i for i in range(n-1,-1, -1)]
        print("qubs = ",qubs)
        self.tt = TernaryTree(self.num_spin_orbitals)

        def append_maj_exc(maj_exc_par, tt, ls: tuple):
            coef = 1
            maj = [tt.get_majorana(i) for i in ls]
            for el in maj:
                coef = coef * el[1]
            maj = [el[0] for el in maj]

            maj = tuple(maj)
            if ls in maj_exc_par:
                prod = maj[0]
                for pauli in maj[1:]:
                    prod, sign = prod_pauli_strings(prod, pauli)
                    coef = coef * sign
                prod = ((qubs.index(gate[0]), gate[1]) for gate in prod)
                cirq.pauli(prod, qubs, coef * maj_exc_par[ls])
                maj_exc_par.pop(ls)

        def single_prep1():
            for i in range(n // 4):
                pauli = self.tt.branch_transposition(2 * i, 0, 2 * i + 1, 0)
                cirq.pauli(pauli, qubs)
                pauli = self.tt.branch_transposition(n // 2 + 2 * i, 0, n // 2 + 2 * i + 1, 0)
                cirq.pauli(pauli, qubs)

        def single_prep2():
            for i in range((n + 2) // 4 - 1):
                pauli = self.tt.branch_transposition(2 * i + 1, 0, 2 * i + 2, 0)
                cirq.pauli(pauli, qubs)
                pauli = self.tt.branch_transposition(n // 2 + 2 * i + 1, 0, n // 2 + 2 * i + 2, 0)
                cirq.pauli(pauli, qubs)

        def single_layer(maj_exc_par, tt):
            for i in range(n):
                ls = sorted([tt[qubs[i]][0].num, tt[qubs[i]][1].num])
                append_maj_exc(maj_exc_par, tt, tuple(ls))

        def double_preparation():
            num = []
            for el in [[2 * i + 1, 2 * i + 2, 2 * i + n + 1, 2 * i + n + 2] for i in range(n // 2)]:
                num = num + el
            # TODO
            for l in range(n // 2 - 1):
                for i in range(l + 1):
                    pauli = self.tt.branch_transposition(n // 2 - l + 2 * i - 1, 0, n // 2 - l + 2 * i, 0)
                    cirq.pauli(pauli, qubs)
                    pauli = self.tt.branch_transposition(n // 2 - l + 2 * i - 1, 1, n // 2 - l + 2 * i, 1)
                    cirq.pauli(pauli, qubs)

            for i in range(self.tt.n_qubits // 2):
                pauli = self.tt.branch_transposition(2 * i, 1, 2 * i + 1, 0)
                cirq.pauli(pauli, qubs)

        def double_antipreparation():
            num = []
            for el in [[2 * i + 1, 2 * i + 2, 2 * i + n + 1, 2 * i + n + 2] for i in range(n // 2)]:
                num = num + el

            for i in reversed(range(self.tt.n_qubits // 2)):
                pauli = self.tt.branch_transposition(2 * i, 1, 2 * i + 1, 0)
                cirq.pauli(pauli, qubs)

            # TODO
            for l in reversed(range(n // 2 - 1)):
                for i in reversed(range(l + 1)):
                    pauli = self.tt.branch_transposition(n // 2 - l + 2 * i - 1, 0, n // 2 - l + 2 * i, 0)
                    cirq.pauli(pauli, qubs=qubs)
                    pauli = self.tt.branch_transposition(n // 2 - l + 2 * i - 1, 1, n // 2 - l + 2 * i, 1)
                    cirq.pauli(pauli,qubs=qubs)

        def swap():
            for i in range(n // 2):
                qubs[2 * i], qubs[2 * i + 1] = qubs[2 * i + 1], qubs[2 * i]
                cirq.swap(2 * i, 2 * i + 1)
            for i in range(n // 2 - 1):
                qubs[2 * i + 1], qubs[2 * i + 2] = qubs[2 * i + 2], qubs[2 * i + 1]
                cirq.swap(2 * i + 1, 2 * i + 2)

        def double_layer(maj_exc_par, tt):
            for i in range(n // 2):
                ls = sorted([tt[qubs[2 * i + 1]][1].num, tt[qubs[2 * i + 1]][0].num,
                             tt[qubs[2 * i]][1].num, tt[qubs[2 * i]][0].num])
                append_maj_exc(maj_exc_par, tt, tuple(ls))

            for i in range(n // 2 - 1):
                ls = sorted([tt[qubs[2 * i + 2]][0].num, tt[qubs[2 * i + 2]][1].num,
                             tt[qubs[2 * i + 1]][0].num, tt[qubs[2 * i + 1]][1].num])
                append_maj_exc(maj_exc_par, tt, tuple(ls))

        N = self.mol.num_spatial_orbitals
        for i in range(n):
            cirq.id(i)
        for i in range(self.num_alpha):
            cirq.x(qubs[i])
        for i in range(self.num_alpha,N - self.num_alpha + 1):
            cirq.id(qubs[i])
        for i in range(self.num_alpha):
            cirq.x(qubs[i + N])
        for i in range(self.num_alpha,N - self.num_alpha + 1):
            cirq.id(qubs[i + N])
        # qubs = list(reversed(qubs))
        print(self.tt)
        for _ in range(n // 2):
            single_layer(self.maj_exc_par, self.tt)
            single_prep1()
            single_prep2()

        double_preparation()

        for _ in range(n // 2):
            double_layer(self.maj_exc_par, self.tt)
            swap()

        for i in range(n // 2):
            pauli = self.tt.branch_transposition(qubs[2 * i], 1, qubs[2 * i + 1], 1)
            cirq.pauli(pauli, qubs)

        for _ in range(n // 2):
            double_layer(self.maj_exc_par, self.tt)
            swap()
        print("maj_exc = ", self.maj_exc_par)
        double_antipreparation()
        return cirq
