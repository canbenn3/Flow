#include "kinematic-wave-model-liu.cpp"
#include <cstdlib>
#include <cstdio>

void timestep_forward(
    GridCell* precomputed_cell_geometry,
    double* water_levels_in,
    double* water_levels_out,
    int width,
    int height
) {
    // Reset all output water levels to zero (we will be adding, not overwriting)
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            water_levels_out[y*width + x] = 0.0;
        }
    }

    // Iterate through each cell and update neighbors
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Discharge discharge = compute_discharge(
                precomputed_cell_geometry[y*width + x],
                water_levels_in[y*width + x]
            );

            int north_y = y - 1;
            if (north_y >= 0) {
                water_levels_out[north_y*width + x] += discharge.north;
            }

            int east_x = x + 1;
            if (east_x < width) {
                water_levels_out[y*width + east_x] += discharge.east;
            }

            int south_y = y + 1;
            if (south_y < height) {
                water_levels_out[south_y*width + x] += discharge.south;
            }

            int west_x = x - 1;
            if (west_x >= 0) {
                water_levels_out[y*width + west_x] += discharge.west;
            }
        }
    }
}

void compute_grid_geometry(
    int* elevations,
    int elevations_width,
    int elevations_height,
    GridCell* geometry
) {
    for (int y = 0; y < elevations_height - 1; y++) {
        for (int x = 0; x < elevations_width - 1; x++) {
            geometry[y*elevations_width + x] = flow_vector_direction(
                elevations[y*elevations_width + x],
                elevations[y*elevations_width + x+1],
                elevations[(y+1)*elevations_width + x+1],
                elevations[(y+1)*elevations_width + x]
            );
        }
    }
}

void run_simulation(
    GridCell* precomputed_cell_geometry,
    double* a,
    double* b,
    int num_timesteps,
    int width,
    int height
) {
    for (int i = 0; i < num_timesteps; i++) {
        if (i % 2 == 0) {
            timestep_forward(precomputed_cell_geometry, a, b, width, height);
        } else {
            timestep_forward(precomputed_cell_geometry, b, a, width, height);
        }
    }
}

int main() {
    int elevations_width = 1024;
    int elevations_height = 1024;
    int* elevations = (int*) malloc(elevations_width*elevations_height * sizeof(int));

    // TODO: read in elevation data from dataset

    // Compute flow steepness and direction for each cell in grid
    int grid_width = elevations_width - 1;
    int grid_height = elevations_height - 1;
    GridCell* geometry = (GridCell*) malloc(grid_width*grid_height * sizeof(GridCell));
    compute_grid_geometry(elevations, elevations_width, elevations_height, geometry);
    free(elevations);

    // Run the simulation
    double* a = (double*) malloc(grid_width*grid_height * sizeof(double));
    double* b = (double*) malloc(grid_width*grid_height * sizeof(double));
    run_simulation(geometry, a, b, 1000, grid_width, grid_height);
    free(b);

    double* result = a;
    
    // TODO: write result out to file
}
