#include "kinematic-wave-model-liu.cpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "read-geo-data.cpp"
#include "write-bmp.cpp"

__global__ void timestep_forward_kernel(
    GridCell *precomputed_cell_geometry,
    double *water_levels_in,
    double *water_levels_out,
    int width,
    int height,
    double dt
)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x >= width || y >= height) return;

    double current_water = water_levels_in[y * width + x];

    // Compute discharge using the new dt parameter
    Discharge discharge = compute_discharge(
        precomputed_cell_geometry[y * width + x],
        current_water,
        dt);

    // Subtract the total water leaving this cell
    double total_outflow = discharge.north + discharge.south + discharge.east + discharge.west;
    water_levels_out[y * width + x] -= total_outflow;

    // Distribute the outflow to the appropriate neighbors
    int north_y = y - 1;
    if (north_y >= 0)
    {
        water_levels_out[north_y * width + x] += discharge.north;
    }

    int east_x = x + 1;
    if (east_x < width)
    {
        water_levels_out[y * width + east_x] += discharge.east;
    }

    int south_y = y + 1;
    if (south_y < height)
    {
        water_levels_out[south_y * width + x] += discharge.south;
    }

    int west_x = x - 1;
    if (west_x >= 0)
    {
        water_levels_out[y * width + west_x] += discharge.west;
    }
}

__global__ void compute_grid_geometry_kernel(
    double *elevations,
    int elevations_width,
    int elevations_height,
    GridCell *geometry
)
{
    int grid_width = elevations_width - 1;
    int grid_height = elevations_height - 1;

    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x < grid_width && y < grid_height)
    {
        geometry[y * grid_width + x] = flow_vector_direction(
            elevations[y * elevations_width + x],
            elevations[y * elevations_width + x + 1],
            elevations[(y + 1) * elevations_width + x + 1],
            elevations[(y + 1) * elevations_width + x]);
    }
}

__global__ void add_rain_kernel(
    double *water_levels,
    int width,
    int height,
    double rain
)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x < grid_width && y < grid_height)
    {
        water_levels[y * width + x] += rain;
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
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if (x < elev_width && y < elev_height)
    {
        int grid_width = elev_width - 1;
        int grid_height = elev_height - 1;

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

        // The new "terrain" the water sees is the actual dirt + the accumulated water
        surface_elevations_out[y * elev_width + x] = terrain_elevations[y * elev_width + x] + avg_water_depth;
    }
}

void run_simulation(
    double *elevation,
    double *surface_elevations,
    double *a,
    double *b,
    int num_timesteps,
    int elev_width,
    int elev_height,
    double dt,
    double total_rainfall_inches
)
{
    int block_size_x = 32;
    int block_size_y = 32;

    int cell_width = elev_width - 1;
    int cell_height = elev_height - 1;

    dim3 blocksPerElevationGrid(ceil(elev_width/(double)block_size_x), ceil(elev_height/(double)block_size_y));
    dim3 blocksPerCellGrid(ceil(cell_width/(double)block_size_x), ceil(cell_height/(double)block_size_y));
    dim3 threadsPerBlock(block_size_x, block_size_y);

    // Add water each timestep to simulate rainfall
    double rain_per_timestep = total_rainfall_inches / num_timesteps;
    GridCell *geometry = (GridCell *)malloc(cell_width * cell_height * sizeof(GridCell));
    for (int i = 0; i < num_timesteps; i++)
    {
        double *in = i % 2 == 0 ? a : b;
        double *out = i % 2 == 0 ? b : a;
        add_rain(in, cell_width, cell_height, rain_per_timestep);
        compute_water_surface_elevations(elevation, in, surface_elevations, elev_width, elev_height);
        compute_grid_geometry(surface_elevations, elev_width, elev_height, geometry);
        timestep_forward(geometry, in, out, cell_width, cell_height, dt);
    }
    free(geometry);
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

    auto [elevations, width, height] = read_geo_data(filename);
    if (elevations == nullptr)
    {
        return 1;
    }

    // Compute flow steepness and direction for each cell in grid
    int grid_width = width - 1;
    int grid_height = height - 1;

    // Setup water grids
    double *a = (double *)calloc(grid_width * grid_height, sizeof(double)); // calloc zeroes the memory
    double *b = (double *)calloc(grid_width * grid_height, sizeof(double));

    double *surface_elevations = (double *)malloc(width * height * sizeof(double));

    // Simulation Parameters
    double dt = 1; // Time step (seconds)

    // Run the simulation
    run_simulation(elevations, surface_elevations, a, b, num_timesteps, width, height, dt, inch_to_meter(rain_inches_total));

    // If num_timesteps is even, 'a' holds the final state. If odd, 'b' holds it.
    double *result = (num_timesteps % 2 == 0) ? a : b;

    write_bmp(filename, grid_width, grid_height, result);

    // Cleanup
    free(a);
    free(b);
    free(surface_elevations);
    free(elevations);

    return 0;
}