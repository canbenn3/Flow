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

__device__ Discharge compute_discharge_for_cell(
    int2 pos,
    double *water_surface_elevations,
    int2 elev_dimens,
    double *water_levels,
    int2 grid_dimens,
    double dt
)
{
    GridCell cell_geometry = flow_vector_direction(
        water_surface_elevations[pos.y * elev_dimens.x + pos.x],
        water_surface_elevations[pos.y * elev_dimens.x + pos.x + 1],
        water_surface_elevations[(pos.y + 1) * elev_dimens.x + pos.x + 1],
        water_surface_elevations[(pos.y + 1) * elev_dimens.x + pos.x]
    );

    double current_water = water_levels[pos.y * grid_dimens.x + pos.x];

    Discharge discharge = compute_discharge(
        cell_geometry,
        current_water,
        dt
    );

    return discharge;
}

__global__ void timestep_forward_kernel(
    double *water_surface_elevations,
    int2 elev_dimens,
    double *water_levels_in,
    double *water_levels_out,
    int2 grid_dimens,
    double dt
)
{
    int2 pos = {
        .x = blockDim.x * blockIdx.x + threadIdx.x,
        .y = blockDim.y * blockIdx.y + threadIdx.y
    };

    if (pos.x < grid_dimens.x && pos.y < grid_dimens.y)
    {
        // Compute inflow by visiting each neighbor
        double inflow = 0;
        for (u_int8_t i = 0; i < 9; i++) {
            int8_t del_x = (i % 3) - 1;
            int8_t del_y = (i / 3) - 1;
            
            // Ignore corners and center
            if ((del_x == 0) ^ (del_y == 0)) {
                int2 neighbor_pos = {pos.x+del_x, pos.y+del_y};
                Discharge neighbor_discharge = compute_discharge_for_cell(
                    neighbor_pos,
                    water_surface_elevations,
                    elev_dimens,
                    water_levels_in,
                    grid_dimens,
                    dt
                );

                // North neighbor: south discharge flows into cell
                if (del_y == -1) {
                    inflow += neighbor_discharge.south;
                }

                // East neighbor: west discharge flows into cell
                if (del_x == 1) {
                    inflow += neighbor_discharge.west;
                }

                // South neighbor: north discharge flows into cell
                if (del_y == 1) {
                    inflow += neighbor_discharge.north;
                }

                // West neighbor: east discharge flows into cell
                if (del_x == -1) {
                    inflow += neighbor_discharge.east;
                }
            }
        }

        // Compute outflow
        Discharge discharge = compute_discharge_for_cell(
            pos,
            water_surface_elevations,
            elev_dimens,
            water_levels_in,
            grid_dimens,
            dt
        );
        double outflow = discharge.north + discharge.south + discharge.east + discharge.west;

        water_levels_out[pos.y * grid_dimens.x + pos.x] += inflow - outflow;
    }
}

__global__ void add_rain_kernel(
    double *water_levels,
    int2 grid_dimens,
    double rain
)
{
    int2 pos = {
        .x = blockDim.x * blockIdx.x + threadIdx.x,
        .y = blockDim.y * blockIdx.y + threadIdx.y
    };

    if (pos.x < grid_dimens.x && pos.y < grid_dimens.y)
    {
        water_levels[pos.y * grid_dimens.x + pos.x] += rain;
    }
}

// A helper function to run at the start of every timestep
__global__ void compute_water_surface_elevations_kernel(
    double *terrain_elevations,
    int2 elev_dimens,
    double *water_levels,
    int2 grid_dimens,
    double *surface_elevations_out
)
{
    int2 pos = {
        .x = blockDim.x * blockIdx.x + threadIdx.x,
        .y = blockDim.y * blockIdx.y + threadIdx.y
    };

    if (pos.x < elev_dimens.x && pos.y < elev_dimens.y)
    {
        double total_water = 0.0;
        int cells_counted = 0;

        // Check the 4 cells touching this corner node.
        // If they are within bounds, add their water depth to our average.

        // Northwest cell
        if (pos.x > 0 && pos.y > 0)
        {
            total_water += water_levels[(pos.y - 1) * grid_dimens.x + (pos.x - 1)];
            cells_counted++;
        }
        // Northeast cell
        if (pos.x < grid_dimens.x && pos.y > 0)
        {
            total_water += water_levels[(pos.y - 1) * grid_dimens.x + pos.x];
            cells_counted++;
        }
        // Southwest cell
        if (pos.x > 0 && pos.y < grid_dimens.y)
        {
            total_water += water_levels[pos.y * grid_dimens.x + (pos.x - 1)];
            cells_counted++;
        }
        // Southeast cell
        if (pos.x < grid_dimens.x && pos.y < grid_dimens.y)
        {
            total_water += water_levels[pos.y * grid_dimens.x + pos.x];
            cells_counted++;
        }

        double avg_water_depth = (cells_counted > 0) ? (total_water / cells_counted) : 0.0;

        // // The new "terrain" the water sees is the actual dirt + the accumulated water
        surface_elevations_out[pos.y * elev_dimens.x + pos.x] = terrain_elevations[pos.y * elev_dimens.x + pos.x] + avg_water_depth;
    }
}

int run_simulation(
    double *elevations_in,
    int2 elev_dimens,
    double *water_levels_out,
    int num_timesteps,
    double dt,
    double total_rainfall_inches
)
{
    int block_size_x = 32;
    int block_size_y = 32;

    int2 grid_dimens = {
        .x = elev_dimens.x - 1,
        .y = elev_dimens.y - 1
    };

    // Add water each timestep to simulate rainfall
    double rain_per_timestep = total_rainfall_inches / num_timesteps;

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
    int grid_width = elev_width - 1;
    int grid_height = elev_height - 1;
    int2 elev_dimens = {elev_width, elev_height};

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
