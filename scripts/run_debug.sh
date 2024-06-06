#!/bin/bash
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
##SBATCH --partition=gp_debug
#SBATCH --qos=gp_debug
#SBATCH --account=bsc08
#SBATCH -t 02:00:00
#SBATCH -o output-%j
#SBATCH -e error-%j
#SBATCH --exclusive


export OMP_DISPLAY_ENV=true
export OMP_NUM_THREADS=48
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

module purge
#module load gcc/8.1.0 openmpi/4.0.1
#module load gcc/11.2.0 openmpi/3.1.1
#module load gcc/11.2.0_binutils openmpi/4.1.2
#module load gcc/7.2.0 openmpi/3.1.1 DDT
#module load gcc/10.1.0 openmpi/3.1.1 DDT/19.0.5
module load gcc/13.2.0 openmpi/4.1.5-gcc ddt

#module load openmp gcc/8.1.0 openmpi/3.1.1 DDT
#openmpi/2.1.3
#openmpi/4.0.1
#module load DDT

make $1
#g++ -O3 -o $1 $1.cpp -fopenmp
#mpiexec -n 2 ./heterogeneity.exe
#ddt --connect mpirun -n 4 ./diff
#ddt --connect
mpirun -n 12  ./tridiagonals 2000 10 1 48 ./test.csv 1> tridiagonals.log 2> tridiagonals.log
