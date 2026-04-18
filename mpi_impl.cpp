/**
 * read-geo-data.cpp requires C++17 or greater.
 * Compile with:  
 * Run with:
 * 
 * MPI implementation of the serial kinematic wave simulation
 * 
 * We use row decomposition becuase we will be using halo rows/ghost cells
 * to so that boudries can communicate. Neighboring cells and rows need to 
 * read from each other to predict flow and everything else.
 * If we used 2D grids, figuring out indexing and exchange would be awful.
 */

#include "kinematic-wave-model-liu.cpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <chrono>
#include <mpi.h>
#include "read-geo-data.cpp"
#include "write-bmp.cpp"
#include "write-jpg.cpp"
#include "util.cpp"

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
    int local_rows,
    int start_row,
    int rank)
{
    // Add water each timestep to simulate rainfall
    double rain_per_timestep = total_rainfall_inches / num_timesteps;
    int grid_width = elev_width - 1;
    GridCell *geometry = (GridCell *)malloc(grid_width * local_rows * sizeof(GridCell));
    for (int i = 0; i < num_timesteps; i++)
    {
        double *in = i % 2 == 0 ? a : b;
        double *out = i % 2 == 0 ? b : a;
        add_rain(in, grid_width, local_rows, rain_per_timestep);
        compute_water_surface_elevations(elevation, in, surface_elevations, elev_width, local_rows+1);
        compute_grid_geometry(surface_elevations, elev_width, local_rows+1, geometry);
        timestep_forward(geometry, in, out, grid_width, local_rows, dt);
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
    std::cout << "USAGE: ./serial <elev.tif> <num_iter> <total_rainfall(inches)>\n";
}

int main(int argc, char *argv[])
{
    //set up and init MPIO  
    int comm_sz, my_rank;
    MPI_Comm comm;
    comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    //have all processes read cmd line for ease - very little data is taken up, probably
    if (argc < 2)
    {
        usage();
        MPI_Finalize();
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

    //set up some global vars
    double *elevations = nullptr;
    double *a = nullptr;
    double *b = nullptr;
    double *surface_elevations = nullptr;
    int width = 0;
    int height = 0;

    //set up termination variable incase we cannot read file
    int read_okay = 1;

    //rank 0 reads the file because it is more expensive if every process reads the file. 
    //Also not every process needs to know elevations
    if (my_rank == 0)
    {
        GeoData data = read_geo_data(filename);
        
        elevations = data.elevations;
        width = data.width;
        height = data.height;
        
        if (elevations == nullptr)
        {
            read_okay = 0;
        }
    }

    MPI_Bcast(&read_okay, 1, MPI_INT, 0, comm);

    if (read_okay == 0)
    {
        MPI_Finalize();
        return 1;
    }

    //broadcast width and height
    MPI_Bcast(&width, 1, MPI_INT, 0, comm);
    MPI_Bcast(&height, 1, MPI_INT, 0, comm);
    
    //set up decomposition math for the water grid 
    int grid_width = width - 1;
    int grid_height = height - 1;

    //find rows per process and remainder rows
    int rows_per_process = grid_height / comm_sz;
    int remainder = grid_height % comm_sz;


    //create sendcounts and displacements for grid
    int *gridSendCounts = (int *)malloc(comm_sz * sizeof(int));
    int *gridDisplacements = (int *)malloc(comm_sz * sizeof(int));

    //create send and displacement for terrian
    int *terrianSendCounts = (int *)malloc(comm_sz * sizeof(int));
    int *terrianDisplacements = (int *)malloc(comm_sz * sizeof(int));
    

    //fill out both grid and terrian send/displacements
    for(int r = 0; r < comm_sz; r++)
    {
        int rows_for_r = rows_per_process + (r < remainder ? 1 : 0);
        int start_row_r = r * rows_per_process + std::min(r, remainder);

        gridSendCounts[r] = rows_for_r * grid_width;
        gridDisplacements[r] = start_row_r * grid_width;

        terrianSendCounts[r] = (rows_for_r + 1) * width;
        terrianDisplacements[r] = start_row_r * width;
    }

    //find grid local_rows, start_row
    int local_rows = gridSendCounts[my_rank] / grid_width; //handles overflow in case grid_height % comm_sz != 0
    int start_row = gridDisplacements[my_rank] / grid_width;

     //find terian local rows because we need more rows and items per row for terrian arrays, cause they are 1 bigger. 
    int terrian_local_rows = local_rows + 1;
    int terrian_start_row = terrianDisplacements[my_rank] / width;


    //calculate if grid rank needs 0, 1 or 2 ghost rows
    //We use ghost_rows because chunks will need to read their neighbors - either the one above or the one below them to calculate movement
    //we wil lhave ghost rows live at the top and bottom of the owned arrays for ease of exchange
    //which means that offsets are used to figure out where to put and pull data from
    int top_ghost_rows = 0;
    int bottom_ghost_rows = 0;

    if (my_rank != 0)
    {
        top_ghost_rows = 1;
    }
    if (my_rank != comm_sz - 1)
    {
        bottom_ghost_rows = 1;
    }

    int ghost_rows = top_ghost_rows + bottom_ghost_rows;
    int owned_grid_offset = top_ghost_rows * grid_width;

    //have root initalize global variables 
    if (my_rank == 0)
    {
        a = (double *)calloc(grid_width * grid_height, sizeof(double));
        b = (double *)calloc(grid_width * grid_height, sizeof(double));
        surface_elevations = (double *)malloc(width * height * sizeof(double));
    }
    

    //set up local grid vars to scatter into
    double *local_a = (double *)calloc(grid_width * (local_rows + ghost_rows), sizeof(double)); // calloc zeroes the memory
    double *local_b = (double *)calloc(grid_width * (local_rows + ghost_rows), sizeof(double));
    
    //set up local terrian vars to scatter into
    double *local_surface_elevations = (double *)malloc(width * terrian_local_rows * sizeof(double));
    double *local_elevations = (double *)malloc(width * terrian_local_rows * sizeof(double));

    //scatter data into a and b. We are ghosting so local_a and local_b are offset by a number of elements to account for the top ghost row
    MPI_Scatterv(a, gridSendCounts, gridDisplacements, MPI_DOUBLE, local_a + owned_grid_offset, gridSendCounts[my_rank], MPI_DOUBLE, 0, comm);
    MPI_Scatterv(b, gridSendCounts, gridDisplacements, MPI_DOUBLE, local_b + owned_grid_offset, gridSendCounts[my_rank], MPI_DOUBLE, 0, comm);

    //scatter data in local_surface 
    MPI_Scatterv(elevations, terrianSendCounts, terrianDisplacements, MPI_DOUBLE, local_elevations, terrianSendCounts[my_rank], MPI_DOUBLE, 0, comm);

    // Simulation Parameters
    double dt = 1; // Time step (seconds)

    // Run the simulation
    auto start_time = std::chrono::high_resolution_clock::now();
    run_simulation(local_elevations, local_surface_elevations, local_a, local_b, num_timesteps, width, terrian_local_rows, dt, inch_to_meter(rain_inches_total), local_rows, start_row, my_rank);
    
    //gather local a and local b into a and b
    double *local_result = (num_timesteps % 2 == 0) ? local_a : local_b;
    double *result = (num_timesteps % 2 == 0) ? a : b ;
    MPI_Gatherv(local_result + owned_grid_offset, gridSendCounts[my_rank], MPI_DOUBLE, result, gridSendCounts, gridDisplacements, MPI_DOUBLE, 0, comm);

    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;
    //if rank zero, print message, and write_jpg
    if (my_rank == 0)
    {
        std::cout << "Simulation completed in " << diff.count() << " seconds." << std::endl;
    

        // write_bmp(filename, grid_width, grid_height, result);
        write_jpg(filename, grid_width, grid_height, result);
    }

    // Cleanup
    free(a);
    free(b);
    free(surface_elevations);
    free(elevations);
    free(local_a);
    free(local_b);
    free(local_surface_elevations);
    free(local_elevations);
    free(gridSendCounts);
    free(gridDisplacements);
    free(terrianSendCounts);
    free(terrianDisplacements);

    MPI_Finalize();

    return 0;
}