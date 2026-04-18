/**
 * read-geo-data.cpp equires C++17 or greater.
 * Compile with `nvcc --std=c++17 cuda.cu -o cuda -lgdal -lstdc++fs`
 */

#include "kernels.cu"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "read-geo-data.cpp"
#include "write-bmp.cpp"

int run_simulation(
    double *elevations_in,
    uint2 elev_dimens,
    double *water_levels_out,
    int num_timesteps,
    double dt,
    double total_rainfall_meters
)
{
    int block_size_x = 32;
    int block_size_y = 32;

    uint2 grid_dimens = {
        .x = elev_dimens.x - 1,
        .y = elev_dimens.y - 1
    };

    // Add water each timestep to simulate rainfall
    double rain_per_timestep = total_rainfall_meters / num_timesteps;

    double *elevations_d;
    double *surface_elevations_d;
    double *water_levels_a_d;
    double *water_levels_b_d;

    cudaError_t ret;
    int elev_len = elev_dimens.x*elev_dimens.y*sizeof(double);
    int grid_len = grid_dimens.x*grid_dimens.y*sizeof(double);

    ret = cudaMalloc((void **) &elevations_d, elev_len);
    if (ret != cudaSuccess) {
        printf("Device memory allocation failed; ret=%d\n", ret);
        return EXIT_FAILURE;
    }

    ret = cudaMemcpy(elevations_d, elevations_in, elev_len, cudaMemcpyHostToDevice);
    if (ret != cudaSuccess) {
        printf("Host to device memory copy failed; ret=%d\n", ret);
        return EXIT_FAILURE;
    }

    ret = cudaMalloc((void **) &surface_elevations_d, elev_len);
    if (ret != cudaSuccess) {
        printf("Device memory allocation failed; ret=%d\n", ret);
        return EXIT_FAILURE;
    }

    ret = cudaMalloc((void **) &water_levels_a_d, grid_len);
    if (ret != cudaSuccess) {
        printf("Device memory allocation failed; ret=%d\n", ret);
        return EXIT_FAILURE;
    }

    ret = cudaMalloc((void **) &water_levels_b_d, grid_len);
    if (ret != cudaSuccess) {
        printf("Device memory allocation failed; ret=%d\n", ret);
        return EXIT_FAILURE;
    }

    dim3 blocksPerElevationGrid(ceil(elev_dimens.x/(double)block_size_x), ceil(elev_dimens.y/(double)block_size_y));
    dim3 blocksPerCellGrid(ceil(grid_dimens.x/(double)block_size_x), ceil(grid_dimens.y/(double)block_size_y));
    dim3 threadsPerBlock(block_size_x, block_size_y);

    for (int i = 0; i < num_timesteps; i++)
    {
        double *in = i % 2 == 0 ? water_levels_a_d : water_levels_b_d;
        double *out = i % 2 == 0 ? water_levels_b_d : water_levels_a_d;

        add_rain_kernel<<<blocksPerCellGrid, threadsPerBlock>>>(
            in, grid_dimens, rain_per_timestep
        );

        ret = cudaDeviceSynchronize();
        if (ret != cudaSuccess) {
            fprintf(stderr, "CUDA synchronize failed; ret=%d\n", ret);
            return EXIT_FAILURE;
        }

        compute_water_surface_elevations_kernel<<<blocksPerElevationGrid, threadsPerBlock>>>(
            elevations_d, elev_dimens, in, grid_dimens, surface_elevations_d
        );

        ret = cudaDeviceSynchronize();
        if (ret != cudaSuccess) {
            fprintf(stderr, "CUDA synchronize failed; ret=%d\n", ret);
            return EXIT_FAILURE;
        }

        timestep_forward_kernel<<<blocksPerCellGrid, threadsPerBlock>>>(
            surface_elevations_d, elev_dimens, in, out, grid_dimens, dt
        );

        ret = cudaDeviceSynchronize();
        if (ret != cudaSuccess) {
            fprintf(stderr, "CUDA synchronize failed; ret=%d\n", ret);
            return EXIT_FAILURE;
        }

        if (i % 10 == 0) {
            fprintf(stderr, "Iteration %d\n", i);
        }
    }

    double *result = (num_timesteps % 2 == 0) ? water_levels_a_d : water_levels_b_d;
    ret = cudaMemcpy(water_levels_out, result, grid_len, cudaMemcpyDeviceToHost);
    if (ret != cudaSuccess) {
        printf("Device to host memory copy failed; ret=%d\n", ret);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

void usage()
{
    std::cout << "USAGE: ./cuda <elev.tif> <num_iter> <total_rainfall(inches)>\n";
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        usage();
        return 1;
    }
    char *filename = argv[1];
    int num_timesteps = 1000; // defaults
    double rain_inches_total = 1;
    if (argc > 2)
    {
        num_timesteps = std::stoi(argv[2]);
    }
    if (argc > 3)
    {
        rain_inches_total = std::stod(argv[3]);
    }

    auto [elevations, elev_width, elev_height] = read_geo_data(filename);
    if (elevations == nullptr)
    {
        return EXIT_FAILURE;
    }
    unsigned int grid_width = (unsigned int) elev_width - 1;
    unsigned int grid_height = (unsigned int) elev_height - 1;
    uint2 elev_dimens = {(unsigned int) elev_width, (unsigned int) elev_height};

    printf("Loaded geo data; %dx%d\n", elev_width, elev_height);

    // Simulation Parameters
    double dt = 1; // Time step (seconds)

    double *water_levels = (double *) malloc(grid_width*grid_height*sizeof(double));

    // Run the simulation
    int ret = run_simulation(
        elevations,
        elev_dimens,
        water_levels,
        num_timesteps,
        dt,
        inch_to_meter(rain_inches_total)
    );
    if (ret != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    write_bmp(filename, grid_width, grid_height, water_levels);

    // Cleanup
    free(water_levels);

    return 0;
}
