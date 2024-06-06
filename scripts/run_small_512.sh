#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --threads-per-core=2
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc08
#SBATCH -t 2-00:00:00
#SBATCH -o small-512-mt_v2.log
#output-%j
#SBATCH -e small-512-mt_v2.log
#error-%j
#SBATCH --exclusive


export OMP_DISPLAY_ENV=true
#export OMP_NUM_THREADS=48
export OMP_PROC_BIND=close
#export OMP_PLACES=threads

module purge
module load gcc/13.2.0 openmpi/4.1.5-gcc

#Small side 2000 substrates 1
for threads in {112..224} 
do
	export OMP_NUM_THREADS=$threads
	#export OMP_PLACES="{0:56},{112:56}" 
	srun --cpus-per-task=112 --threads-per-core=2 ./tridiagonals_512 5000 10 2 1 ./apply_v2/vectormt_v2/c_${cpus}_th_${threads}_512.csv #1> small.log 2> small.log
done 
#Medium 5000 subs 2

#Large 10000 subs 4

#Extra large side 14000 subs 8

#Massive side 20000 subs 8
