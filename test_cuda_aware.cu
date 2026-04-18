/**
 * Compile with `nvcc -c test_cuda_aware.cu -o cuda-test.o && mpicxx cuda-test.o -o cuda-test -lcudart`
 * Tested successfully with cuda/13.1.0, openmpi/5.0.8-gpu, and OMPI_MCA_opal_cuda_support=1
 */

#include <iostream>
#include <mpi.h>
#include <cuda_runtime.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int test_value = (rank == 0) ? 42 : 0;
    int *d_ptr;

    // 1. Allocate on GPU
    cudaError_t err = cudaMalloc(&d_ptr, sizeof(int));
    if (err != cudaSuccess) {
        fprintf(stderr, "Rank %d: cudaMalloc failed\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // 2. Perform CUDA-Aware MPI call
    // We are passing a GPU pointer (d_ptr) directly to MPI
    MPI_Request req;
    MPI_Isend(&test_value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &req);
    MPI_Request_free(&req);
    MPI_Recv(d_ptr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // 3. Verify by copying back to a local host variable
    int result = 0;
    cudaMemcpy(&result, d_ptr, sizeof(int), cudaMemcpyDeviceToHost);

    if (result == 42) {
        printf("Rank %d: Success! Received %d on GPU.\n", rank, result);
    } else {
        printf("Rank %d: Failed! Received %d.\n", rank, result);
    }

    cudaFree(d_ptr);
    MPI_Finalize();
    return 0;
}
