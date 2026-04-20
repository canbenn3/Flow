# Single GPU

The single GPU implementation uses a pull approach to update the water level at each cell by calculating and incrementing by the inflow from neighboring cells. The obvious alternative is a push approach that calculates the outflow from self and then updates each neighbor accordingly. Although the former involves more overhead (flow direction vectors must be computed independently for each neighbor of any given cell or pre-computed and stored in memory), it prevents data race conditions that could occur in the latter approach when two cells that share a neighbor try to write to that neighbor simultaneously.

Elevation map data is read from a file on the host and then copied to the GPU. Water levels and elevation data are stored in global device memory.

### Instructions

To run on CHPC, request an allocation with a GPU. We tested on `kingspeak-gpu` nodes with Nvidia GTX Titan X (SP, 12GB) with CUDA 11.8.0.

1. Request an allocation using the following (or similar) constraints.
```
salloc -N 1 -n 1 --ntasks-per-node=1 --gres=gpu:1 --time=00:15:00 -p kingspeak-gpu -A kingspeak-gpu --gres=gpu
```
2. Load required modules.
```
module load cuda gdal
```
3. Compile.
```
nvcc --std=c++17 cuda.cu -o cuda -lgdal -lstdc++fs
```
4. Run.
```
./cuda <elev.tif> <num_iter> <total_rainfall(inches)>
```

# Distributed Memory GPU

The distributed memory GPU implementation reuses kernels from the regular GPU implementation, following the "pull" strategy described above. Orchestration of GPU nodes is accomplished using OpenMPI.

Elevation map data is read from a file on the root host and then distributed among worker nodes.

To simplify tiling, the elevation map is split along the vertical axis only. Since each cell has a dependency on it's direct neighbors and neighbors-once-removed, a halo region of 2 rows is used. Adjacent GPUs exchange the most up-to-date water levels in the halo region at each time step.

Data is exchanged directly between GPUs, bypassing the host via a CUDA-aware MPI environment with Remote Direct Memory Access (RDMA).

### Instructions

To run on CHPC, request an allocation with more than GPU. Since our implementation simply uses the default GPU known to CUDA, each GPU must live on a distinct node. We tested on `kingspeak-gpu` nodes with Nvidia GTX Titan X (SP, 12GB) with CUDA 11.8.0 and on `notchpeak-gpu` nodes with CUDA 11.8.0.

1. Request an allocation using the following (or similar) constraints.
```
salloc -N 2 -n 2 --ntasks-per-node=1 --gres=gpu:1 --time=00:15:00 -p kingspeak-gpu -A kingspeak-gpu --gres=gpu
```
2. Load required modules.
```
module load cuda/11.8.0 openmpi/4.1.6-gpu gdal
```
3. Compile object files with `nvcc` and link with `mpicxx`.
```
nvcc -c --std=c++17 cuda-distributed.cu -o cuda-distributed.o && mpicxx cuda-distributed.o -o cuda-distributed -lcudart -lgdal -lstdc++fs
```
4. Enable CUDA RDMA.
```
export OMPI_MCA_opal_cuda_support=1
```
5. Remove page lock limits that otherwise prevent initialization of the network adapters used by RDMA.
```
ulimit -l unlimited
```
6. Run.
```
mpiexec -n 2 cuda-distributed <elev.tif> <num_iter> <total_rainfall(inches)>
```
