#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --threads-per-core=2
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc08
#SBATCH -t 1:00:00
#SBATCH -o output-%j
#SBATCH -e error-%j
#SBATCH --exclusive


source load_modules_MN5.sh

export OMP_DISPLAY_ENV=false
export OMP_NUM_THREADS=112
export OMP_PROC_BIND=close
export OMP_PLACES=threads

#Parameters definition
SIDE=5000 #Cube side in um
SUBS=2    #Number of substrates
DELTA=10  #Delta {x,y,z} in um
BLOCKS=1  #Number of blocks to parallelize  in multi nodes
CSV="./test" #Path to store timing 

module load ddt

srun ./$1 $SIDE $DELTA $SUBS $BLOCKS $CSV

#ddt --connect ./correctness $SIDE $DELTA $SUBS $BLOCKS $CSV