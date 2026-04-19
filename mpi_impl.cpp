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
    double *water_levels_in, //includes ghost rows,
    double *water_levels_out, //owned region only,
    int width,
    int height, //of local owned slice - does not count ghost rows,
    int top_ghost_row_count, 
    int bottom_ghost_row_count,
    double dt // timestep
)
{
    //figure out how many total rows
    int total_in_rows = height + top_ghost_row_count + bottom_ghost_row_count;
    
    // Iterate through each cell, subtract outflow, and add inflow to neighbors
    for (int y = 0; y < height; y++)
    {
        int in_y = y + top_ghost_row_count;

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
                    int n_y = in_y + del_y;
                    if (n_x < width && n_x >= 0 && n_y >= 0 && n_y < total_in_rows)
                    {
                        double current_water = water_levels_in[(in_y + del_y) * width + (x + del_x)];

                        Discharge neighbor_discharge = compute_discharge(
                            precomputed_cell_geometry[(in_y + del_y) * width + (x + del_x)],
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
                precomputed_cell_geometry[in_y * width + x],
                water_levels_in[in_y * width + x],
                dt);
            // Subtract the total water leaving this cell
            double outflow = own_discharge.north + own_discharge.south + own_discharge.east + own_discharge.west;
            water_levels_out[y * width + x] = water_levels_in[in_y * width + x] + inflow - outflow;
        }
    }
}

//geometry has ghost rows - we don't need to change anything because we only pass through owned slice
void compute_grid_geometry(
    double *elevations,
    int elevations_width,
    int elevations_height,
    GridCell *geometry) //local owned slice only - techinally has ghost row, but height won't read it)
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
//for mpi, y is the top of the local terrian slice, not the global evelvations
//because we only pass in owned, ghost rows should be safe
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

//helper function to exchange ghost rows between processes. Essentially, this is here to make it so functions can read from each other. Assume it takes full slice, not the owned slice
void halo_exchange_water 
(
    double *water_levels,
    int grid_width,
    int local_rows,
    int rank,
    int comm_sz,
    int owned_grid_offset
)
{   
    //calculate pointers for various parts of the water_levels given local rows and such. 
    double *top_ghost_row = water_levels;
    double *owned_start_row = water_levels + owned_grid_offset;
    
    double *owned_end_row = water_levels + owned_grid_offset + (local_rows - 1) * grid_width;
    double *bottom_ghost_row = water_levels + owned_grid_offset + local_rows * grid_width;

    //determine neighbors
    int top_neighbor = rank - 1;
    int bottom_neighbor = rank + 1;

    //exchange logic
    if(rank != 0)
    {
        //exchange top ghost row with top neighbor. Using MPI_Sendrecv to avoid deadlocks
        MPI_Sendrecv(
            owned_start_row,
            grid_width,
            MPI_DOUBLE,
            top_neighbor,
            0,
            top_ghost_row,
            grid_width,
            MPI_DOUBLE,
            top_neighbor,
            0,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
        );
    }

    if(rank != comm_sz - 1)
    {
        //exchange bottom ghost row with bottom neighbor. using MPI_Sendrecv to avoid deadlocks
        MPI_Sendrecv(
            owned_end_row,
            grid_width,
            MPI_DOUBLE,
            bottom_neighbor,
            0,
            bottom_ghost_row,
            grid_width,
            MPI_DOUBLE,
            bottom_neighbor,
            0,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
        );
    }

}
//helper function to exchange rows of gridcell grids 
//- is needed instead of water cause structs are different data type than doubles
//we are sending bytes instead.
void halo_exchange_geometry(
    GridCell *geometry, //full local slice, not just owned slice,
    int grid_width,
    int local_rows,
    int rank,
    int comm_sz, 
    int owned_grid_offset
)
{
    //calculate pointers for various parts of the geometery given local rows and such. 
    GridCell *top_ghost_row = geometry;
    GridCell *owned_start_row = geometry + owned_grid_offset;
    
    GridCell *owned_end_row = geometry + owned_grid_offset + (local_rows - 1) * grid_width;
    GridCell *bottom_ghost_row = geometry + owned_grid_offset + local_rows * grid_width;

    //determine neighbors
    int top_neighbor = rank - 1;
    int bottom_neighbor = rank + 1;

    //how many bytes
    int row_byte_count = grid_width * sizeof(GridCell);

    //if valid, send and recive top rows
    if(rank != 0)
    {
        MPI_Sendrecv(
            owned_start_row, 
            row_byte_count,
            MPI_BYTE,
            top_neighbor,
            0,
            top_ghost_row,
            row_byte_count,
            MPI_BYTE,
            top_neighbor,
            0,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
        );
    }

    //if valid, send and receive bottom rows
    if(rank != comm_sz - 1)
    {
        MPI_Sendrecv(
            owned_end_row, 
            row_byte_count,
            MPI_BYTE,
            bottom_neighbor,
            0,
            bottom_ghost_row,
            row_byte_count,
            MPI_BYTE,
            bottom_neighbor,
            0,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
        );
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
    int rank,
    int p)
{
    //compute things necessary for the simulation
    double rain_per_timestep = total_rainfall_inches / num_timesteps;
    int grid_width = elev_width - 1;
    int top_ghost_row_count = 0;
    int bottom_ghost_row_count = 0;
    //find ghost rows
    if (rank != 0) 
    {
        top_ghost_row_count = 1;
    }
    if (rank != p - 1)
    {
        bottom_ghost_row_count = 1;
    }

    //find grid offset for ghost rows
    int owned_grid_offset = top_ghost_row_count * grid_width;
    //run simulation.
    GridCell *geometry = (GridCell *)malloc(grid_width * (local_rows + top_ghost_row_count + bottom_ghost_row_count) * sizeof(GridCell));
    for (int i = 0; i < num_timesteps; i++)
    {
        double *in = i % 2 == 0 ? a : b;
        double *out = i % 2 == 0 ? b : a;

        //use offset to find owned region start:
        double *owned_in = in + owned_grid_offset;
        double *owned_out = out + owned_grid_offset;
        GridCell *owned_geometry = geometry + owned_grid_offset;

        add_rain(owned_in, grid_width, local_rows, rain_per_timestep);
        halo_exchange_water(in, grid_width, local_rows, rank, p, owned_grid_offset);
        compute_water_surface_elevations(elevation, owned_in, surface_elevations, elev_width, elev_height);

        compute_grid_geometry(surface_elevations, elev_width, elev_height, owned_geometry);

        //once we know the grid geometry, we must echange top and bottom rows of the owned slices
        halo_exchange_geometry(geometry, grid_width, local_rows, rank, p, owned_grid_offset);
        timestep_forward(geometry, in, owned_out, grid_width, local_rows, top_ghost_row_count, bottom_ghost_row_count, dt);
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

    //find grid local_rows
    int local_rows = gridSendCounts[my_rank] / grid_width; //handles overflow in case grid_height % comm_sz != 0;

     //find terian local rows because we need more rows and items per row for terrian arrays, cause they are 1 bigger. 
    int terrian_local_rows = local_rows + 1;


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

    //have root initalize global variables - we don't need globabl surface_elevations or goemtery because those are mostly just done on locakl
    if (my_rank == 0)
    {
        a = (double *)calloc(grid_width * grid_height, sizeof(double));
        b = (double *)calloc(grid_width * grid_height, sizeof(double));
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
    run_simulation(local_elevations, local_surface_elevations, local_a, local_b, num_timesteps, width, terrian_local_rows, dt, inch_to_meter(rain_inches_total), local_rows, my_rank, comm_sz);
    
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