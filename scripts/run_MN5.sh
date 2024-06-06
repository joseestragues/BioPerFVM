#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc08
#SBATCH -t 2-00:00:00
#SBATCH -o small-256.log
#output-%j
#SBATCH -e small-256.log
#error-%j
#SBATCH --exclusive


export OMP_DISPLAY_ENV=false
#export OMP_NUM_THREADS=48
export OMP_PROC_BIND=close
export OMP_PLACES=threads

module purge
module load gcc/13.2.0 openmpi/4.1.5-gcc

#Small side 2000 substrates 1
for threads in {1..112}
do
	export OMP_NUM_THREADS=$threads
	srun --cpus-per-task=$threads ./tridiagonals 5000 10 2 1 ./apply_v2/vector/${threads}_threads.csv #1> small.log 2> small.log
done 
#Medium 5000 subs 2

#Large 10000 subs 4

#Extra large side 14000 subs 8

#Massive side 20000 subs 8
