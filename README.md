# Wannier90 Hamiltonian Tools
Simple tools used to manipulate Wannier90 Hamiltonians written by Alexandru B. Georgescu. Please cite the work below if you use it; please note that these tools are under constant improvement and are not meant as production-level code at this point. 

The simple matlab scripts used here should make it quite easy to verify the accuracy of a Wannier Hamiltonian in fitting a DFT calculation, as well as perform various types of post-processing. This includes, for example, rotation of the orbital basis to obtain a locally diagonal density matrix, which can be particularly important in certain materials, for example in dihalides and trihalides, as shown in our work here: (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.245153)

Example 1 for TiCl2 is used to illustrate the diagonalization and rotation procedure discussed in our work here: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.245153 . 

The tools provided can be used as follows. First we want to get the bands out of the bands.out file, as obtained from pw.x, in a format we can load in Matlab. We can use:

for python 2:
getbands.py bands.out > bands.dat
for python3
getbands3.py bands.out > bands.dat

The script plotbandsfromthis will load the bands.dat file, and the hr file from the Wannier90 output. This can then be used to plot the DFT band structure, the Wannier bands as well as Wannier orbital projected bands. PDOSrotated can then be used to plot the projected density of states of the orbitals, either rotated, or not. The comments inside the scripts should make the usage self-explanatory. In order to obtain the rotated Wannier basis (the .vesta files), one simply linearly mixes the .xsf files obtained from the Wannierization using Vesta's isosurface tools.



This code was written by Alexandru B. Georgescu, while working in collaboration with James Rondinelli and Andrew Millis, and builds on previous code written with Sohrab Ismail-Beigi used to manipulate Wannier functions as an input to slave-boson methods as described here: 'Boson Subsidiary Solver (BoSS) v1.1', AB Georgescu, M Kim, S Ismail-Beigi, Computer Physics Communications 265, 107991 ; the getbands.py script was written by Andrei Malashevich.
