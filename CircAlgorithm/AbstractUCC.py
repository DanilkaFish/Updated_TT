from __future__ import annotations

import openfermion
from qiskit.circuit import Parameter, QuantumCircuit
from openfermion import transforms
from openfermionpyscf import run_pyscf
from qiskit_nature.second_q.drivers import PySCFDriver

from .utils import simplification, single_exc, double_exc
from abc import ABC, abstractmethod
from qiskit_nature.second_q.problems import ElectronicBasis, ElectronicStructureProblem


class AbstractUCC(ABC):
    def __init__(self,
                 geometry: str= 'H 0 0 0; H 0 0 0.7414',
                 basis: str = "sto-3g",
                 multiplicity: int = 1,
                 charge: int = 0,
                 description: str = '',
                 filename: str = '',
                 data_directory: str = None
                 ):
        driver = PySCFDriver(atom=geometry,
                              basis=basis,
                              spin= (multiplicity - 1)//2,
                              charge=charge)
        self.mol = driver.run()
        self.fermionic_op = self.mol.hamiltonian.second_q_op()
        self.n_qubits = int(self.mol.num_spin_orbitals)
        self.maj_exc_par = self.get_excitations()

    @property
    def num_alpha(self):
        return self.mol.num_alpha

    @property
    def num_beta(self):
        return self.mol.num_beta

    @property
    def num_spatial_orbitals(self):
        return self.n_qubits//2

    @property
    def num_spin_orbitals(self):
        return self.n_qubits

    def get_excitations(self) -> {(int, int) | (int, int)}:
        """
        Method to get parametrized UpUCCSDG excitations via via majorana operators
        """
        ladder_exc_par = {}

        ladder_exc = self.get_alpha_excitations()
        ladder_exc += self.get_beta_excitations()
        ladder_exc += self.get_double_excitations()

        for op in ladder_exc:
            ladder_exc_par[op] = Parameter('t_' + ','.join([str(i) for i in op]))
            

        maj_exc_par = {}
        for op in ladder_exc_par:
            if len(op) == 2:
                for new_op, sign in single_exc(op):
                    maj_exc_par[new_op] = sign * ladder_exc_par[op] + maj_exc_par.get(new_op, 0)
            if len(op) == 4:
                for new_op, sign in double_exc(op):
                    maj_exc_par[new_op] = sign * ladder_exc_par[op] + maj_exc_par.get(new_op, 0)

        maj_exc_par = simplification(maj_exc_par)
        return maj_exc_par

    @abstractmethod
    def get_alpha_excitations(self) -> list[tuple[int, int]]:
        pass

    @abstractmethod
    def get_beta_excitations(self) -> list[tuple[int, int]]:
        pass

    @abstractmethod
    def get_double_excitations(self) -> list[tuple[int, int, int, int]]:
        pass

    @abstractmethod
    def get_parametrized_circuit(self) -> QuantumCircuit:
        pass

