#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc08
#SBATCH -t 12:00:00
#SBATCH -o output-%j
#SBATCH -e error-%j
#SBATCH --exclusive


export OMP_DISPLAY_ENV=true
#export OMP_NUM_THREADS=48
export OMP_PROC_BIND=close
export OMP_PLACES=threads

module purge
module load gcc/13.2.0 openmpi/4.1.5-gcc

#Small side 2000 substrates 1
for threads in 32 64 96 112
do
	
	echo "MPI run execution with ${threads}"
	export OMP_NUM_THREADS=$threads
	mpirun -n 1 --map-by ppr:1:node ./tridiagonals 5000 10 2 1 ./outputs/mpirun_${threads}_threads.csv
	echo "SRUN execution with ${threads}"
	export OMP_NUM_THREADS=$threads
	srun --cpus-per-task=$threads ./tridiagonals 5000 10 2 1 ./outputs/srun_${threads}_threads.csv #1> small.log 2> small.log
done 
#Medium 5000 subs 2

#Large 10000 subs 4

#Extra large side 14000 subs 8

#Massive side 20000 subs 8
