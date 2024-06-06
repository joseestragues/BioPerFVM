#!/bin/bash
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --qos=bsc_ls
#SBATCH -t 1-00:00:00
#SBATCH -o test-16
#SBATCH -e /dev/null
#SBATCH --exclusive


export OMP_DISPLAY_ENV=false
export OMP_NUM_THREADS=48
export OMP_PROC_BIND=close
export OMP_PLACES=threads

module purge
#module load openmp gcc/7.2.0 openmpi/3.1.1 # gcc/8.1.0 openmpi/3.1.1  #openmpi/2.1.3  #openmpi/3.1.4 #openmpi/2.1.3
#openmpi/4.0.1
#module load DDT/18.2
module load gcc/8.1.0 openmpi/3.1.1 DDT/19.0.5


export nodes=16

for voxels in 1000 1500 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000
do
	for substrates in  1 2 4 8 12 16 24 32
	do
		for factor in 1 2 4
		do
			overlap=$((nodes * factor))
			halfside=$((voxels * 2))
			echo "Executing HF: ${voxels} | Substrates: ${substrates} | Steps: ${overlap}"
			srun ./tridiagonals $halfside 10 $substrates $overlap ./timing/voxels_${voxels}/substrates_${substrates}/factor_${factor}/${nodes}_node.csv
		done
	done	
done
