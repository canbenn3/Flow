#include "kinematic-wave-model-liu.cpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "read-geo-data.cpp"
#include "write-bmp.cpp"

void timestep_forward(
    GridCell *precomputed_cell_geometry,
    double *water_levels_in,
    double *water_levels_out,
    int width,
    int height,
    double dt // timestep
)
{
    // Initialize output array with the current water levels.
    // We must keep the water that *doesn't* move, then add/subtract the changes.
    std::memcpy(water_levels_out, water_levels_in, width * height * sizeof(double));

    // Iterate through each cell, subtract outflow, and add inflow to neighbors
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            double current_water = water_levels_in[y * width + x];

            // Optimization: Skip calculations if the cell is essentially dry
            if (current_water < 1e-6)
                continue;

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
    }
}

void compute_grid_geometry(
    double *elevations,
    int elevations_width,
    int elevations_height,
    GridCell *geometry)
{
    int grid_width = elevations_width - 1;

    for (int y = 0; y < elevations_height - 1; y++)
    {
        for (int x = 0; x < elevations_width - 1; x++)
        {
            geometry[y * grid_width + x] = flow_vector_direction(
                elevations[y * elevations_width + x],
                elevations[y * elevations_width + x + 1],
                elevations[(y + 1) * elevations_width + x + 1],
                elevations[(y + 1) * elevations_width + x]);
        }
    }
}

void add_rain(double *water_levels, int width, int height, double rain)
{
    for (int i = 0; i < width * height; i++)
    {
        water_levels[i] += rain;
    }
}

void run_simulation(
    GridCell *precomputed_cell_geometry,
    double *a,
    double *b,
    int num_timesteps,
    int width,
    int height,
    double dt, // Added time step
    double total_rainfall_inches)
{
    // Add water each timestep to simulate rainfall
    double rain_per_timestep = total_rainfall_inches / num_timesteps;
    for (int i = 0; i < num_timesteps; i++)
    {
        if (i % 2 == 0)
        {
            add_rain(a, width, height, rain_per_timestep);
            timestep_forward(precomputed_cell_geometry, a, b, width, height, dt);
        }
        else
        {
            add_rain(b, width, height, rain_per_timestep);
            timestep_forward(precomputed_cell_geometry, b, a, width, height, dt);
        }
    }
}

void printGrid(double arr[], int size, int cols)
{
    for (int i = 0; i < size; ++i)
    {
        // Print element with width 4
        std::cout << std::setw(4) << arr[i] << " ";

        // Print newline after every 'cols' elements
        if ((i + 1) % cols == 0)
        {
            std::cout << std::endl;
        }
    }
    std::cout << std::endl; // Final newline
}

void usage()
{
    std::cout << "USAGE: ./serial <elev.tif> <num_iter> <total_rainfall(inches)>\n";
}

double inch_to_meter(double inch)
{
    return inch * .0254;
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
        int num_timesteps = std::stoi(argv[2]);
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

    printGrid(elevations, width * height, width);

    // Compute flow steepness and direction for each cell in grid
    int grid_width = width - 1;
    int grid_height = height - 1;
    GridCell *geometry = (GridCell *)malloc(grid_width * grid_height * sizeof(GridCell));
    compute_grid_geometry(elevations, width, height, geometry);
    free(elevations); // We don't need elevations anymore; geometry holds the static gradients!

    // Setup water grids
    double *a = (double *)calloc(grid_width * grid_height, sizeof(double)); // calloc zeroes the memory
    double *b = (double *)calloc(grid_width * grid_height, sizeof(double));

    // Simulation Parameters
    double dt = 1; // Time step (seconds)

    // Run the simulation
    run_simulation(geometry, a, b, num_timesteps, grid_width, grid_height, dt, inch_to_meter(rain_inches_total));

    // If num_timesteps is even, 'a' holds the final state. If odd, 'b' holds it.
    double *result = (num_timesteps % 2 == 0) ? a : b;

    std::cout << "final water level\n";
    printGrid(result, grid_width * grid_height, grid_width);
    write_bmp(filename, grid_width, grid_height, result);

    // Cleanup
    free(a);
    free(b);
    free(geometry);

    return 0;
}