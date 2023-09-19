# deterministic_exact_unitary_inversion

This repository contains

- Qiskit code of a quantum circuit universally reversing qubit-unitary operation deterministically and exactly
- SDP code to obtain the optimal success probability or fidelity of unitary inversion

These codes are accompanied to the following paper:

- Satoshi Yoshida, Akihito Soeda, and Mio Murao, Reversing unknown qubit-unitary operation, deterministically and exactly, [Phys. Rev. Lett. 131, 120602 (2023)](https://doi.org/10.1103/PhysRevLett.131.120602), [arXiv:2209.02907](https://arxiv.org/abs/2209.02907).

## Requirement

The [Qiskit](https://qiskit.org) code of the qubit-unitary inversion is written in jupyter notebook and requires the packages [qiskit](https://doi.org/10.5281/zenodo.2573505) and [pylatexenc](https://pypi.org/project/pylatexenc).

The SDP code to obtain the optimal unitary inversion is written in Matlab and requires the following interpreter:

- [CVX](http://cvxr.com): a Matlab-based convex modeling framework

These codes also use functions of QETLAB ([QETLAB](https://qetlab.com): A MATLAB Toolbox for Quantum Entanglement), but all used functions are contained in the subfolder [QETLAB_used_functions](https://github.com/sy3104/deterministic_exact_unitary_inversion/tree/main/QETLAB_used_functions).

It has been tested on MATLAB R2021b and CVX 2.2.

## Description

The main scripts of this repository are

- [qubit_unitary_inversion_circuit.ipynb](https://github.com/sy3104/deterministic_exact_unitary_inversion/blob/main/qubit_unitary_inversion_circuit.ipynb): Qiskit code to generate a quantum circuit of a deterministic exact qubit-unitary inversion for a random qubit-unitary operation and a random input qubit state.  In the quantum circuit shown in this page, the qubits are rearranged from the quantum circuit shown in the paper to reduce the number of SWAP gates.  However, the quantum circuit shown in this page is equivalent to that shown in the paper.
- [run_optimal_untiary_inversion.m](https://github.com/sy3104/deterministic_exact_unitary_inversion/blob/main/run_optimal_untiary_inversion.m): Code to obtain the maximal maximal fidelity of transforming n calls of any 'd'-dimensional unitary operations into its inverse map.

The script [run_optimal_untiary_inversion.m](https://github.com/sy3104/deterministic_exact_unitary_inversion/blob/main/run_optimal_untiary_inversion.m) makes use of the functions in the subfolder [QETLAB_used_functions](https://github.com/sy3104/deterministic_exact_unitary_inversion/tree/main/QETLAB_used_functions) and the following Matlab codes in the subfolder [sdp](https://github.com/sy3104/deterministic_exact_unitary_inversion/tree/main/sdp):

- [deterministic_parallel_unitary_inversion.m](https://github.com/sy3104/deterministic_exact_unitary_inversion/blob/main/sdp/deterministic_parallel_unitary_inversion.m): Code to obtain the maximal fidelity of transforming parallel n calls of any d-dimensional unitary operations into its inverse map.
- [deterministic_sequential_unitary_inversion.m](https://github.com/sy3104/deterministic_exact_unitary_inversion/blob/main/sdp/deterministic_sequential_unitary_inversion.m): Code to obtain the maximal fidelity of transforming sequential n calls of any d-dimensional unitary operations into its inverse map.

These Matlab codes utilize group-theoretic calculations done by a [Sagemath](https://www.sagemath.org) code shown in [young_diagrams.sage](https://github.com/sy3104/deterministic_exact_unitary_inversion/blob/main/sdp/young_diagrams.sage).


## License

This code is under the [MIT license](https://opensource.org/licenses/MIT).
