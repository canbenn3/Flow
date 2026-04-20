#!/bin/bash
#SBATCH --time=00:30:00          # Walltime (HH:MM:SS) - adjust as needed
#SBATCH --nodes=1                # Number of nodes
#SBATCH --ntasks-per-node=1      # Number of tasks (keep at 1 for OpenMP)
#SBATCH --cpus-per-task=1       # Number of threads (OpenMP cores)
#SBATCH --partition=kingspeak    # Partition name
#SBATCH --account=usucs6030    # Replace with your actual account name
#SBATCH --job-name=flow_sim
#SBATCH --output=logs/%j.out     # Standard output log (%j is job ID)
#SBATCH --error=logs/%j.err      # error log
#SBATCH --mem=32G

module load gcc/8.5.0            
module load gdal

mkdir -p output
mkdir -p logs

# USAGE: ./mp <thread_count> <elev.tif> <num_iter> <total_rainfall(inches)>
./mp 16 ../elevation-maps/cache_valley.tif 100 25
./mp 16 ../elevation-maps/USU.tif 100 25
./mp 16 ../elevation-maps/Hyrum.tif 100 25
./mp 16 ../elevation-maps/Hospital.tif 100 25
