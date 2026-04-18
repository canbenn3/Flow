/**
 * read-geo-data.cpp equires C++17 or greater.
 * Compile with `nvcc --std=c++17 cuda.cu -o cuda -lgdal -lstdc++fs`
 */

#include "kinematic-wave-model-liu.cpp"
#include <cstdlib>

__device__ Discharge compute_discharge_for_cell(
    uint2 pos,
    double *water_surface_elevations,
    uint2 elev_dimens,
    double *water_levels,
    uint2 grid_dimens,
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
    uint2 elev_dimens,
    double *water_levels_in,
    double *water_levels_out,
    uint2 grid_dimens,
    double dt
)
{
    uint2 pos = {
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
                uint2 neighbor_pos = {pos.x+del_x, pos.y+del_y};
                if (neighbor_pos.x < grid_dimens.x && neighbor_pos.y < grid_dimens.y)
                {
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
    uint2 grid_dimens,
    double rain
)
{
    uint2 pos = {
        .x = blockDim.x * blockIdx.x + threadIdx.x,
        .y = blockDim.y * blockIdx.y + threadIdx.y
    };

    if (pos.x < grid_dimens.x && pos.y < grid_dimens.y)
    {
        water_levels[pos.y * grid_dimens.x + pos.x] += rain;
    }
}

__global__ void compute_water_surface_elevations_kernel(
    double *terrain_elevations,
    uint2 elev_dimens,
    double *water_levels,
    uint2 grid_dimens,
    double *surface_elevations_out
)
{
    uint2 pos = {
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
