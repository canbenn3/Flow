#!/bin/sh

#SBATCH --time=00-00:15:00

#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH -p kingspeak-gpu
#SBATCH -A kingspeak-gpu

#SBATCH -o slurmjob-%j@%N.out #stdout file in format slurmjob-SLURM_JOB_ID.out-NODEID

cd /uufs/chpc.utah.edu/common/home/u6074057/classes/cs5030/Flow

echo "---INIT ENVIRONMENT---"
module load cuda/11.8.0 openmpi/4.1.6-gpu gdal && export OMPI_MCA_opal_cuda_support=1 && ulimit -l unlimited

echo "---COMPILE---"
nvcc -c --std=c++17 cuda-distributed.cu -o cuda-distributed.o && mpicxx cuda-distributed.o -o cuda-distributed -lcudart -lgdal -lstdc++fs

echo "---Hospital.tif 100 25 (2 GPUs)---"
time mpiexec -n 2 cuda-distributed ../elevation-maps/Hospital.tif 100 25

echo "---USU.tif 100 25 (2 GPUs)---"
time mpiexec -n 2 cuda-distributed ../elevation-maps/USU.tif 100 25

echo "---Hyrum.tif 100 25 (2 GPUs)---"
time mpiexec -n 2 cuda-distributed ../elevation-maps/Hyrum.tif 100 25

echo "---cache_valley.tif 100 25 (2 GPUs)---"
time mpiexec -n 2 cuda-distributed ../elevation-maps/cache_valley.tif 100 25
