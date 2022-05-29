# Wannier90HamiltonianTools
Simple tools used to manipulate Wannier90 Hamiltonians written by Alexandru B. Georgescu. Please cite the work below if you use it; please note that these tools are under constant improvement and are not meant as production-level code at this point. 

The simple matlab scripts used here should make it quite easy to verify the accuracy of a Wannier Hamiltonian in fitting a DFT calculation, as well as perform various types of post-processing. This includes, for example, rotation of the orbital basis to obtain a locally diagonal density matrix, which can be particularly important in certain materials, for example in dihalides and trihalides, as shown in our work here: https://arxiv.org/abs/2110.04665

Example 1 is used to illustrate the diagonalization and rotation procedure discussed in our work here: https://arxiv.org/abs/2110.04665 . The script plotbandsfromthis will plot the DFT band structure, the Wannier bands as well as Wannier 'fat-bands' projected to certain orbitals. PDOSrotated can then be used to plot the projected density of states of the orbitals, either rotated, or not. The comments inside the scripts should make the usage self-explanatory.



This code was written by Alexandru B. Georgescu, while working in collaboration with James Rondinelli and Andrew Millis, and builds on previous code written with Sohrab Ismail-Beigi used to manipulate Wannier functions as an input to slave-boson methods as described here: 'Boson Subsidiary Solver (BoSS) v1.1', AB Georgescu, M Kim, S Ismail-Beigi, Computer Physics Communications 265, 107991 ; the getbands.py script was written by Andrei Malashevich.
