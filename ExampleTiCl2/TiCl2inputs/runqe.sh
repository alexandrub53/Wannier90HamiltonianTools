#!/bin/bash                                                                                                                                                           
#SBATCH -J SIO                                                                                                                                                  
#SBATCH --ntasks 40                                                                                                                                                   
#SBATCH --partition=bnl                                                                                                                                               
#SBATCH -o out.%J                                                                                                                                                     
#SBATCH -e err.%J                                                                                                                                                     
#SBATCH -t 72:00:00                                                                                                                                                   
#SBATCH --mail-type=ALL                                                                                                                                               
module purge
module add gcc openmpi/1.10.7-hfi   lib/openblas/0.2.19-haswell   lib/fftw3/3.3.6-pl1 slurm

#mpirun -np 40  /mnt/home/ageorgescu/qe/qe-6.6/bin/pw.x -nk 40 -ndiag 1 < scf.in > scf.out
#mpirun -np 40  /mnt/home/ageorgescu/qe/qe-6.6/bin/pw.x -nk 40 -ndiag 1 < bands.in > bands.out
#mpirun -np 40  /mnt/home/ageorgescu/qe/qe-6.6/bin/pw.x -nk 40 -ndiag 1 < nscf.in > nscf.out
mpirun -np 40 /mnt/home/ageorgescu/qe/qe-6.6/bin/wannier90.x --pp nno
mpirun -np 40  /mnt/home/ageorgescu/qe/qe-6.6/bin/pw2wannier90.x < nno.pw2wan.in > nno.pw2wan.out
mpirun -np 28 /mnt/home/ageorgescu/qe/qe-6.6/bin/wannier90.x  nno
