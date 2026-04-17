/**
 * read-geo-data.cpp equires C++17 or greater.
 * Compile with `nvcc --std=c++17 cuda.cu -o cuda -lgdal -lstdc++fs`
 */

#include "kinematic-wave-model-liu.cpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "read-geo-data.cpp"
#include "write-bmp.cpp"

__global__ void timestep_forward_kernel(
    double *water_surface_elevations,
    double *water_levels_in,
    double *water_levels_out,
    int grid_width,
    int grid_height,
    double dt
)
{
    int elevations_width = grid_width + 1;

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x < grid_width && y < grid_height)
    {
        GridCell cell_geometry = flow_vector_direction(
            water_surface_elevations[y * elevations_width + x],
            water_surface_elevations[y * elevations_width + x + 1],
            water_surface_elevations[(y + 1) * elevations_width + x + 1],
            water_surface_elevations[(y + 1) * elevations_width + x]
        );

        double current_water = water_levels_in[y * grid_width + x];

        Discharge discharge = compute_discharge(
            cell_geometry,
            current_water,
            dt
        );

        // Subtract the total water leaving this cell
        double total_outflow = discharge.north + discharge.south + discharge.east + discharge.west;
        water_levels_out[y * grid_width + x] -= total_outflow;

        // Distribute the outflow to the appropriate neighbors
        int north_y = y - 1;
        if (north_y >= 0)
        {
            water_levels_out[north_y * grid_width + x] += discharge.north;
        }

        int east_x = x + 1;
        if (east_x < grid_width)
        {
            water_levels_out[y * grid_width + east_x] += discharge.east;
        }

        int south_y = y + 1;
        if (south_y < grid_height)
        {
            water_levels_out[south_y * grid_width + x] += discharge.south;
        }

        int west_x = x - 1;
        if (west_x >= 0)
        {
            water_levels_out[y * grid_width + west_x] += discharge.west;
        }
    }
}

__global__ void add_rain_kernel(
    double *water_levels,
    int grid_width,
    int grid_height,
    double rain
)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x < grid_width && y < grid_height)
    {
        water_levels[y * grid_width + x] += rain;
    }
}

// A helper function to run at the start of every timestep
__global__ void compute_water_surface_elevations_kernel(
    double *terrain_elevations,
    double *water_levels,
    double *surface_elevations_out,
    int elev_width,
    int elev_height
)
{
    int grid_width = elev_width - 1;
    int grid_height = elev_height - 1;

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x < elev_width && y < elev_height)
    {
        double total_water = 0.0;
        int cells_counted = 0;

        // Check the 4 cells touching this corner node.
        // If they are within bounds, add their water depth to our average.

        // Northwest cell
        if (x > 0 && y > 0)
        {
            total_water += water_levels[(y - 1) * grid_width + (x - 1)];
            cells_counted++;
        }
        // Northeast cell
        if (x < grid_width && y > 0)
        {
            total_water += water_levels[(y - 1) * grid_width + x];
            cells_counted++;
        }
        // Southwest cell
        if (x > 0 && y < grid_height)
        {
            total_water += water_levels[y * grid_width + (x - 1)];
            cells_counted++;
        }
        // Southeast cell
        if (x < grid_width && y < grid_height)
        {
            total_water += water_levels[y * grid_width + x];
            cells_counted++;
        }

        double avg_water_depth = (cells_counted > 0) ? (total_water / cells_counted) : 0.0;

        // // The new "terrain" the water sees is the actual dirt + the accumulated water
        surface_elevations_out[y * elev_width + x] = terrain_elevations[y * elev_width + x] + avg_water_depth;
    }
}

int run_simulation(
    double *elevations_in,
    int elev_width,
    int elev_height,
    double *water_levels_out,
    int num_timesteps,
    double dt,
    double total_rainfall_inches
)
{
    int block_size_x = 32;
    int block_size_y = 32;

    int grid_width = elev_width - 1;
    int grid_height = elev_height - 1;

    // Add water each timestep to simulate rainfall
    double rain_per_timestep = total_rainfall_inches / num_timesteps;

    double *elevations_d;
    double *surface_elevations_d;
    double *water_levels_a_d;
    double *water_levels_b_d;

    cudaError_t ret;
    int elev_len = elev_width*elev_height*sizeof(double);
    int grid_len = grid_width*grid_height*sizeof(double);

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

    dim3 blocksPerElevationGrid(ceil(elev_width/(double)block_size_x), ceil(elev_height/(double)block_size_y));
    dim3 blocksPerCellGrid(ceil(grid_width/(double)block_size_x), ceil(grid_height/(double)block_size_y));
    dim3 threadsPerBlock(block_size_x, block_size_y);

    for (int i = 0; i < num_timesteps; i++)
    {
        double *in = i % 2 == 0 ? water_levels_a_d : water_levels_b_d;
        double *out = i % 2 == 0 ? water_levels_b_d : water_levels_a_d;

        add_rain_kernel<<<blocksPerCellGrid, threadsPerBlock>>>(in, grid_width, grid_height, rain_per_timestep);

        ret = cudaDeviceSynchronize();
        if (ret != cudaSuccess) {
            fprintf(stderr, "CUDA synchronize failed; ret=%d\n", ret);
            return EXIT_FAILURE;
        }

        compute_water_surface_elevations_kernel<<<blocksPerElevationGrid, threadsPerBlock>>>(elevations_d, in, surface_elevations_d, elev_width, elev_height);

        ret = cudaDeviceSynchronize();
        if (ret != cudaSuccess) {
            fprintf(stderr, "CUDA synchronize failed; ret=%d\n", ret);
            return EXIT_FAILURE;
        }

        timestep_forward_kernel<<<blocksPerCellGrid, threadsPerBlock>>>(surface_elevations_d, in, out, grid_width, grid_height, dt);

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
    int grid_width = elev_width - 1;
    int grid_height = elev_height - 1;

    printf("Loaded geo data; %dx%d\n", elev_width, elev_height);

    // Simulation Parameters
    double dt = 1; // Time step (seconds)

    double *water_levels = (double *) malloc(grid_width*grid_height*sizeof(double));

    // Run the simulation
    int ret = run_simulation(elevations, elev_width, elev_height, water_levels, num_timesteps, dt, inch_to_meter(rain_inches_total));
    if (ret != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    write_bmp(filename, grid_width, grid_height, water_levels);

    // Cleanup
    free(water_levels);

    return 0;
}
