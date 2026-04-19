/**
 * read-geo-data.cpp equires C++17 or greater.
 * Configure environment with `module load cuda/11.8.0 openmpi/4.1.6-gpu gdal && export OMPI_MCA_opal_cuda_support=1 && ulimit -l unlimited`
 * Compile with `nvcc -c --std=c++17 cuda-distributed.cu -o cuda-distributed.o && mpicxx cuda-distributed.o -o cuda-distributed -lcudart -lgdal -lstdc++fs`
 */

#include "kernels.cu"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "read-geo-data.cpp"
#include "write-bmp.cpp"
#include <mpi.h>

// #define USE_HOST_COPYBACK_BUF

void usage()
{
    std::cout << "USAGE: ./cuda-distributed <elev.tif> <num_iter> <total_rainfall(inches)>\n";
}

int main(int argc, char *argv[])
{
    int comm_size;
    int rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char* filename;
    int num_timesteps;
    double total_rainfall_meters;
    double* master_elevations = NULL;
    unsigned int master_grid_width, master_grid_height, process_grid_height;
    int offsets[comm_size];
    int grid_rows[comm_size];

    if (rank == 0)
    {
        if (argc < 2)
        {
            usage();
            MPI_Finalize();
            return 0;
        }
        filename = argv[1];
        num_timesteps = 1000; // defaults
        double rain_inches_total = 1;
        if (argc > 2)
        {
            num_timesteps = std::stoi(argv[2]);
        }
        if (argc > 3)
        {
            rain_inches_total = std::stod(argv[3]);
        }
        total_rainfall_meters = inch_to_meter(rain_inches_total);


        auto [elevations, elev_width, elev_height] = read_geo_data(filename);
        if (elevations == NULL)
        {
            return EXIT_FAILURE;
        }
        master_elevations = elevations;

        master_grid_width = (unsigned int) elev_width - 1;
        master_grid_height = (unsigned int) elev_height - 1;

        printf("Loaded geo data; %dx%d\n", elev_width, elev_height);

        unsigned int grid_rows_per_process = master_grid_height / comm_size;

        int offset_cursor = 0;
        for (int i = 0; i < (comm_size - 1); i++) {
            offsets[i] = offset_cursor;
            grid_rows[i] = grid_rows_per_process;
            offset_cursor += grid_rows_per_process;
        }

        // Last process takes the remaining data
        offsets[comm_size - 1] = offset_cursor;
        grid_rows[comm_size - 1] = master_grid_height - offset_cursor;

        // Adjust boundaries so that they overlap (halos)
        for (int i = 0; i < (comm_size - 1); i++) {
            grid_rows[i] += 1;
            offsets[i + 1] -= 1;
            grid_rows[i + 1] += 1;
        }
    }

    MPI_Bcast(&num_timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&total_rainfall_meters, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&master_grid_width, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&master_grid_height, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(grid_rows, 1, MPI_INT, &process_grid_height, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Elevation map is divided along the y-axis
    // Grid width is the same for all processes; height is different
    uint2 grid_dimens_with_halo = {master_grid_width, process_grid_height};
    uint2 elev_dimens_with_halo = {master_grid_width + 1, process_grid_height + 1};

    uint2 grid_dimens = grid_dimens_with_halo;

    int rank_above = rank - 1;
    int rank_below = rank + 1;
    int padding_above = 0;

    if (rank_above >= 0) {
        grid_dimens.y -= 1;
        padding_above = 1;
    }

    if (rank_below < comm_size) {
        grid_dimens.y -= 1;
    }

    // Allocate GPU memory

    double *elevations_d;
    double *surface_elevations_d;
    double *water_levels_a_d;
    double *water_levels_b_d;

    cudaError_t ret;
    int grid_len = grid_dimens.x*grid_dimens.y*sizeof(double);
    int elev_len_with_halo = elev_dimens_with_halo.x*elev_dimens_with_halo.y*sizeof(double);
    int grid_len_with_halo = grid_dimens_with_halo.x*grid_dimens_with_halo.y*sizeof(double);
    int master_grid_len = master_grid_width*master_grid_height*sizeof(double);

    ret = cudaMalloc((void **) &elevations_d, elev_len_with_halo);
    if (ret != cudaSuccess) {
        fprintf(stderr, "Device memory allocation failed; ret=%d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    ret = cudaMalloc((void **) &surface_elevations_d, elev_len_with_halo);
    if (ret != cudaSuccess) {
        fprintf(stderr, "Device memory allocation failed; ret=%d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    ret = cudaMalloc((void **) &water_levels_a_d, grid_len_with_halo);
    if (ret != cudaSuccess) {
        fprintf(stderr, "Device memory allocation failed; ret=%d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    ret = cudaMemset(water_levels_a_d, 0, grid_len_with_halo);
    if (ret != cudaSuccess) {
        fprintf(stderr, "Failed to zero memory; ret=%d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    ret = cudaMalloc((void **) &water_levels_b_d, grid_len_with_halo);
    if (ret != cudaSuccess) {
        fprintf(stderr, "Device memory allocation failed; ret=%d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    // Distribute sections of the elevation map to each node

    if (rank == 0) {
        for (int i = 0; i < comm_size; i++) {
            int elev_offset = offsets[i] * elev_dimens_with_halo.x;
            int elev_len = (grid_rows[i] + 1) * elev_dimens_with_halo.x * sizeof(double);
            MPI_Request req;
            MPI_Isend(master_elevations + elev_offset, elev_len, MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
        }
    }
    MPI_Recv(elevations_d, elev_len_with_halo, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Rank %d: received elevation map data\n", rank);

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        free(master_elevations);
    }

    // Run the simulation on the GPU

    int block_size_x = 32;
    int block_size_y = 32;
    double dt = 1;

    dim3 blocksPerElevationGrid(ceil(elev_dimens_with_halo.x/(double)block_size_x), ceil(elev_dimens_with_halo.y/(double)block_size_y));
    dim3 blocksPerCellGrid(ceil(grid_dimens_with_halo.x/(double)block_size_x), ceil(grid_dimens_with_halo.y/(double)block_size_y));
    dim3 threadsPerBlock(block_size_x, block_size_y);

    double rain_per_timestep = total_rainfall_meters / num_timesteps;

    for (int i = 0; i < num_timesteps; i++)
    {
        double *in = i % 2 == 0 ? water_levels_a_d : water_levels_b_d;
        double *out = i % 2 == 0 ? water_levels_b_d : water_levels_a_d;

        if (rank_above >= 0) {
            MPI_Request req;
            MPI_Isend(in + grid_dimens_with_halo.x, grid_dimens_with_halo.x, MPI_DOUBLE, rank_above, 0, MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (rank_below < comm_size) {
            MPI_Request req;
            MPI_Isend(in + (grid_dimens_with_halo.y - 2) * grid_dimens_with_halo.x, grid_dimens_with_halo.x, MPI_DOUBLE, rank_below, 0, MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
            MPI_Recv(in + (grid_dimens_with_halo.y - 1) * grid_dimens_with_halo.x, grid_dimens_with_halo.x, MPI_DOUBLE, rank_below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (rank_above >= 0) {
            MPI_Recv(in, grid_dimens_with_halo.x, MPI_DOUBLE, rank_above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        add_rain_kernel<<<blocksPerCellGrid, threadsPerBlock>>>(
            in, grid_dimens_with_halo, rain_per_timestep
        );

        ret = cudaDeviceSynchronize();
        if (ret != cudaSuccess) {
            fprintf(stderr, "CUDA synchronize failed; ret=%d\n", ret);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            return EXIT_FAILURE;
        }

        compute_water_surface_elevations_kernel<<<blocksPerElevationGrid, threadsPerBlock>>>(
            elevations_d, elev_dimens_with_halo, in, grid_dimens_with_halo, surface_elevations_d
        );

        ret = cudaDeviceSynchronize();
        if (ret != cudaSuccess) {
            fprintf(stderr, "CUDA synchronize failed; ret=%d\n", ret);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            return EXIT_FAILURE;
        }

        timestep_forward_kernel<<<blocksPerCellGrid, threadsPerBlock>>>(
            surface_elevations_d, elev_dimens_with_halo, in, out, grid_dimens_with_halo, dt
        );

        ret = cudaDeviceSynchronize();
        if (ret != cudaSuccess) {
            fprintf(stderr, "CUDA synchronize failed; ret=%d\n", ret);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            return EXIT_FAILURE;
        }

        if (rank == 0 && i % 10 == 0) {
            printf("Iteration %d\n", i);
        }
    }

    double *result = (num_timesteps % 2 == 0) ? water_levels_a_d : water_levels_b_d;
    double *water_levels = NULL;
    if (rank == 0) {
        water_levels = (double *) malloc(master_grid_len);
        if (water_levels == NULL) {
            fprintf(stderr, "Unable to allocate result buffer on host; len=%d\n", master_grid_len);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            return EXIT_FAILURE;
        }

        // When we scattered data, the boundaries were set to overlap in order
        // to configure halos. When we gather data, we want to copy only the
        // regions without the halos. Remove the overlap.
        for (int i = 0; i < (comm_size - 1); i++) {
            grid_rows[i] -= 1;
            offsets[i + 1] += 1;
            grid_rows[i + 1] -= 1;
        }

        for (int i = 0; i < comm_size; i++) {
            grid_rows[i] *= grid_dimens.x * sizeof(double);
            offsets[i] *= grid_dimens.x * sizeof(double);
        }
    }

    double *result_offset_from_halo = result + (padding_above * grid_dimens_with_halo.x);

    #ifndef USE_HOST_COPYBACK_BUF
    MPI_Gatherv(result_offset_from_halo, grid_len, MPI_BYTE, water_levels, grid_rows, offsets, MPI_BYTE, 0, MPI_COMM_WORLD);
    #endif

    #ifdef USE_HOST_COPYBACK_BUF

    double *water_levels_partial = (double *) malloc(grid_len);
    if (water_levels_partial == NULL) {
        fprintf(stderr, "Unable to allocate copyback buffer on host; len=%d\n", grid_len);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    ret = cudaMemcpy(water_levels_partial, result_offset_from_halo, grid_len, cudaMemcpyDeviceToHost);
    if (ret != cudaSuccess) {
        fprintf(stderr, "Device to host memory copy failed; ret=%d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    MPI_Gatherv(water_levels_partial, grid_len, MPI_BYTE, water_levels, grid_rows, offsets, MPI_BYTE, 0, MPI_COMM_WORLD);
    free(water_levels_partial);

    #endif

    if (rank == 0) {
        write_bmp(filename, master_grid_width, master_grid_height, water_levels);
        free(water_levels);
    }

    MPI_Finalize();

    return 0;
}
