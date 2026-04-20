#include "kinematic-wave-model-liu.cpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <omp.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
// #include "read-geo-data.cpp"
// #include "write-bmp.cpp"
// #include "write-jpg.cpp"
// #include "util.cpp"

// Forward Declarations: Tell the compiler these exist somewhere else
struct GeoData
{
    double *elevations;
    int width;
    int height;
};
GeoData read_geo_data(const char *filename);
double inch_to_meter(double inch);
int write_jpg(const char *filename, int width, int height, double *data);
int write_bmp(const char *filename, int width, int height, double *data);

void timestep_forward(
    GridCell *precomputed_cell_geometry,
    double *water_levels_in,
    double *water_levels_out,
    int width,
    int height,
    double dt // timestep
)
{
    // Iterate through each cell, subtract outflow, and add inflow to neighbors
#pragma omp for
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            // Calculate water flowing into cell (outflow could cause race conditions)
            double inflow = 0;

            for (u_int8_t i = 0; i < 9; i++)
            {
                int8_t del_x = (i % 3) - 1;
                int8_t del_y = (i / 3) - 1;

                // Ignore corners and center
                if ((del_x == 0) ^ (del_y == 0))
                {
                    int n_x = x + del_x;
                    int n_y = y + del_y;
                    if (n_x < width && n_x >= 0 && n_y >= 0 && n_y < height)
                    {
                        double current_water = water_levels_in[(y + del_y) * width + (x + del_x)];

                        Discharge neighbor_discharge = compute_discharge(
                            precomputed_cell_geometry[(y + del_y) * width + (x + del_x)],
                            current_water,
                            dt);

                        // North neighbor: south discharge flows into cell
                        if (del_y == -1)
                        {
                            inflow += neighbor_discharge.south;
                        }

                        // East neighbor: west discharge flows into cell
                        if (del_x == 1)
                        {
                            inflow += neighbor_discharge.west;
                        }

                        // South neighbor: north discharge flows into cell
                        if (del_y == 1)
                        {
                            inflow += neighbor_discharge.north;
                        }

                        // West neighbor: east discharge flows into cell
                        if (del_x == -1)
                        {
                            inflow += neighbor_discharge.east;
                        }
                    }
                }
            }
            Discharge own_discharge = compute_discharge(
                precomputed_cell_geometry[y * width + x],
                water_levels_in[y * width + x],
                dt);
            // Subtract the total water leaving this cell
            double outflow = own_discharge.north + own_discharge.south + own_discharge.east + own_discharge.west;
            water_levels_out[y * width + x] = water_levels_in[y * width + x] + inflow - outflow;
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

#pragma omp for
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
#pragma omp for
    for (int i = 0; i < width * height; i++)
    {
        water_levels[i] += rain;
    }
}

// A helper function to run at the start of every timestep
void compute_water_surface_elevations(
    double *terrain_elevations,
    double *water_levels,
    double *surface_elevations_out,
    int elev_width,
    int elev_height)
{
    int grid_width = elev_width - 1;
    int grid_height = elev_height - 1;

#pragma omp for
    for (int y = 0; y < elev_height; y++)
    {
        for (int x = 0; x < elev_width; x++)
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

            // The new "terrain" the water sees is the actual dirt + the accumulated water
            surface_elevations_out[y * elev_width + x] = terrain_elevations[y * elev_width + x] + avg_water_depth;
        }
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
    double dt, // Added time step
    double total_rainfall_inches,
    int thread_count)
{
    double rain_per_timestep = total_rainfall_inches / num_timesteps;
    int grid_width = elev_width - 1;
    int grid_height = elev_height - 1;
    GridCell *geometry = (GridCell *)malloc(grid_width * grid_height * sizeof(GridCell));
#pragma omp parallel num_threads(thread_count)
    {
        for (int i = 0; i < num_timesteps; i++)
        {
            double *in = i % 2 == 0 ? a : b;
            double *out = i % 2 == 0 ? b : a;
            add_rain(in, grid_width, grid_height, rain_per_timestep);
            compute_water_surface_elevations(elevation, in, surface_elevations, elev_width, elev_height);
            compute_grid_geometry(surface_elevations, elev_width, elev_height, geometry);
            timestep_forward(geometry, in, out, grid_width, grid_height, dt);
        }
    }
    free(geometry);
}

// Helper function for debugging (only use for small .tif files)
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
    std::cout << "USAGE: ./serial <thread_count> <elev.tif> <num_iter> <total_rainfall(inches)>\n";
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        usage();
        return 1;
    }
    int thread_count = std::stoi(argv[1]);
    char *filename = argv[2];
    int num_timesteps = 1000; // defaults
    double rain_inches_total = 1;

    if (argc > 3)
    {
        num_timesteps = std::stoi(argv[3]);
    }
    if (argc > 4)
    {
        rain_inches_total = std::stod(argv[4]);
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
    auto start_time = std::chrono::high_resolution_clock::now();
    run_simulation(elevations, surface_elevations, a, b, num_timesteps, width, height, dt, inch_to_meter(rain_inches_total), thread_count);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;
    std::cout << filename << " completed in " << diff.count() << " seconds with " << thread_count << " threads." << std::endl;

    // If num_timesteps is even, 'a' holds the final state. If odd, 'b' holds it.
    double *result = (num_timesteps % 2 == 0) ? a : b;

    //write_bmp(filename, grid_width, grid_height, result);
    write_jpg(filename, grid_width, grid_height, result);
    // Cleanup
    free(a);
    free(b);
    free(surface_elevations);
    free(elevations);

    return 0;
}