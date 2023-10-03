import numpy as np
from qiskit.quantum_info.operators import Pauli
from .BaseTree import st_enumeration
from .TernaryTree import TernaryTree



def _pauli_table_TT(**kwargs):
    """
    This function is used by InitstateTTInfo for vacuum preparation.
    """
    tt =  kwargs.get("tt", None)
    if tt == None:
        tt = TernaryTree(nmodes)

    n_qubits = tt.n_qubits
    nmodes = tt.nmodes
    branches = tt.branches()
    num_list = tt.enum_list
    print(num_list)
    def branch_to_str(branch):
        pauli_list = ["I"]*n_qubits
        pauli_str = ""
        for el in branch:
            pauli_list[el[0]] = el[1]
        for pauli in pauli_list:
            pauli_str += pauli
        return pauli_str
    pauli_table_str = []
    for i in range(nmodes):
        pauli_table_str.append( (branch_to_str(branches[num_list.index(2*i )]) ,  branch_to_str(branches[num_list.index(2*i + 1)]) ))
    return [pauli_table_str,[[1,-1j] for _ in range(len(pauli_table_str))]]

def pauli_table_TT(**kwargs):
    """
    This function is used by TernaryTreeMapper to obtain pauli_table
    """
    pauli_table = []
    pauli_table_str = _pauli_table_TT(**kwargs)[0]
    for i in range(len(pauli_table_str)):
        pauli_table.append( (Pauli( pauli_table_str[i][0] ), Pauli( pauli_table_str[i][1])) )
    return pauli_table

def _pauli_table_JW(nmodes, num_list = None):
    
    pt = []
    for i in range(nmodes):
        a_z = np.asarray([1] * i + [0] + [0] * (nmodes - i - 1), dtype=bool)
        a_x = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        b_z = np.asarray([1] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        b_x = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        pt.append(["Z" * i + "X" + "I"*(nmodes - i - 1), "Z" * i + "Y" + "I"*(nmodes - i - 1)])
    return [pt,[1,-1j]*len(pt)]

def pauli_table_JW(nmodes,**kwargs):
    pauli_table = []
    for i in range(nmodes):
        a_z = np.asarray([1] * i + [0] + [0] * (nmodes - i - 1), dtype=bool)
        a_x = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        b_z = np.asarray([1] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        b_x = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        # c_z = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        # c_x = np.asarray([0] * nmodes, dtype=bool)
        pauli_table.append((Pauli((a_z, a_x)), Pauli((b_z, b_x))))
    #     print(pauli_table)
    return pauli_table


def pauli_table_BK(nmodes, **kwargs):

    pauli_table = []
    def parity_set(j, n):
        """
        Computes the parity set of the j-th orbital in n modes.

        Args:
            j (int) : the orbital index
            n (int) : the total number of modes

        Returns:
            numpy.ndarray: Array of mode indices
        """
        indices = np.array([])
        if n % 2 != 0:
            return indices

        if j < n / 2:
            indices = np.append(indices, parity_set(j, n / 2))
        else:
            indices = np.append(
                indices, np.append(parity_set(j - n / 2, n / 2) + n / 2, n / 2 - 1)
            )
        return indices

    def update_set(j, n):
        """
        Computes the update set of the j-th orbital in n modes.

        Args:
            j (int) : the orbital index
            n (int) : the total number of modes

        Returns:
            numpy.ndarray: Array of mode indices
        """
        indices = np.array([])
        if n % 2 != 0:
            return indices
        if j < n / 2:
            indices = np.append(indices, np.append(n - 1, update_set(j, n / 2)))
        else:
            indices = np.append(indices, update_set(j - n / 2, n / 2) + n / 2)
        return indices

    def flip_set(j, n):
        """
        Computes the flip set of the j-th orbital in n modes.

        Args:
            j (int) : the orbital index
            n (int) : the total number of modes

        Returns:
            numpy.ndarray: Array of mode indices
        """
        indices = np.array([])
        if n % 2 != 0:
            return indices
        if j < n / 2:
            indices = np.append(indices, flip_set(j, n / 2))
        elif n / 2 <= j < n - 1:
            indices = np.append(indices, flip_set(j - n / 2, n / 2) + n / 2)
        else:
            indices = np.append(
                np.append(indices, flip_set(j - n / 2, n / 2) + n / 2), n / 2 - 1
            )
        return indices

    pauli_table = []
    # FIND BINARY SUPERSET SIZE
    bin_sup = 1
    while nmodes > np.power(2, bin_sup):
        bin_sup += 1
    # DEFINE INDEX SETS FOR EVERY FERMIONIC MODE
    update_sets = []
    update_pauli = []

    parity_sets = []
    parity_pauli = []

    flip_sets = []

    remainder_sets = []
    remainder_pauli = []
    for j in range(nmodes):

        update_sets.append(update_set(j, np.power(2, bin_sup)))
        update_sets[j] = update_sets[j][update_sets[j] < nmodes]

        parity_sets.append(parity_set(j, np.power(2, bin_sup)))
        parity_sets[j] = parity_sets[j][parity_sets[j] < nmodes]

        flip_sets.append(flip_set(j, np.power(2, bin_sup)))
        flip_sets[j] = flip_sets[j][flip_sets[j] < nmodes]

        remainder_sets.append(np.setdiff1d(parity_sets[j], flip_sets[j]))

        update_pauli.append(Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool))))
        parity_pauli.append(Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool))))
        remainder_pauli.append(
            Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool)))
        )
        for k in range(nmodes):
            if np.in1d(k, update_sets[j]):
                update_pauli[j].x[k] = True
            if np.in1d(k, parity_sets[j]):
                parity_pauli[j].z[k] = True
            if np.in1d(k, remainder_sets[j]):
                remainder_pauli[j].z[k] = True

        x_j = Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool)))
        x_j.x[j] = True
        y_j = Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool)))
        y_j.z[j] = True
        y_j.x[j] = True
        pauli_table.append(
            (
                parity_pauli[j] & x_j & update_pauli[j],
                remainder_pauli[j] & y_j & update_pauli[j],
            )
        )

    # PauliList has the phase information unlike deprecated PauliTable.
    # Here, phase is unnecessary, so the following removes phase.
    for pauli1, pauli2 in pauli_table:
        pauli1.phase = 0
        pauli2.phase = 0
    return pauli_table