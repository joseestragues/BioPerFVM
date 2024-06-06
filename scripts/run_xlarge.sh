#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc08
#SBATCH -t 1-00:00:00
#SBATCH -o xlarge.log
#output-%j
#SBATCH -e xlarge.log
#error-%j
#SBATCH --exclusive


export OMP_DISPLAY_ENV=false
#export OMP_NUM_THREADS=48
export OMP_PROC_BIND=close
export OMP_PLACES=threads

module purge
module load gcc/13.2.0 openmpi/4.1.5-gcc

#Small side 2000 substrates 1

#Medium 5000 subs 2

#Large 10000 subs 4

#Extra large side 14000 subs 8
for threads in {1..224}
do
	export OMP_NUM_THREADS=$threads
	srun --cpus-per-task=$threads ./tridiagonals 14000 10 8 1 ./timing_mn5/xl/${threads}_threads.csv #1> small.log 2> small.log
done 
#Massive side 20000 subs 8
