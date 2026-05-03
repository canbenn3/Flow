#include "kinematic-wave-model-liu.cpp"
#include <cstdlib>

__device__ float sample_surface_elevation_node(
    uint2 pos,
    const float *terrain_elevations,
    uint2 elev_dimens,
    const float *water_levels,
    uint2 grid_dimens
)
{
    float total_water = 0.0f;
    int cells_counted = 0;

    if (pos.x > 0 && pos.y > 0)
    {
        total_water += water_levels[(pos.y - 1) * grid_dimens.x + (pos.x - 1)];
        cells_counted++;
    }
    if (pos.x < grid_dimens.x && pos.y > 0)
    {
        total_water += water_levels[(pos.y - 1) * grid_dimens.x + pos.x];
        cells_counted++;
    }
    if (pos.x > 0 && pos.y < grid_dimens.y)
    {
        total_water += water_levels[pos.y * grid_dimens.x + (pos.x - 1)];
        cells_counted++;
    }
    if (pos.x < grid_dimens.x && pos.y < grid_dimens.y)
    {
        total_water += water_levels[pos.y * grid_dimens.x + pos.x];
        cells_counted++;
    }

    const float avg_water_depth = (cells_counted > 0) ? (total_water / cells_counted) : 0.0f;
    return terrain_elevations[pos.y * elev_dimens.x + pos.x] + avg_water_depth;
}

__device__ Discharge compute_discharge_for_cell(
    uint2 pos,
    const float *terrain_elevations,
    uint2 elev_dimens,
    const float *water_levels,
    uint2 grid_dimens,
    double dt
)
{
    GridCell cell_geometry = flow_vector_direction(
        sample_surface_elevation_node({pos.x, pos.y}, terrain_elevations, elev_dimens, water_levels, grid_dimens),
        sample_surface_elevation_node({pos.x + 1, pos.y}, terrain_elevations, elev_dimens, water_levels, grid_dimens),
        sample_surface_elevation_node({pos.x + 1, pos.y + 1}, terrain_elevations, elev_dimens, water_levels, grid_dimens),
        sample_surface_elevation_node({pos.x, pos.y + 1}, terrain_elevations, elev_dimens, water_levels, grid_dimens)
    );

    double current_water = static_cast<double>(water_levels[pos.y * grid_dimens.x + pos.x]);

    Discharge discharge = compute_discharge(
        cell_geometry,
        current_water,
        dt
    );

    return discharge;
}

__global__ void timestep_forward_kernel(
    const float *terrain_elevations,
    uint2 elev_dimens,
    const float *water_levels_in,
    float *water_levels_out,
    uint2 grid_dimens,
    double dt,
    double rain_per_timestep
)
{
    uint2 pos = {
        .x = blockDim.x * blockIdx.x + threadIdx.x,
        .y = blockDim.y * blockIdx.y + threadIdx.y
    };

    if (pos.x < grid_dimens.x && pos.y < grid_dimens.y)
    {
        // Compute inflow by visiting each neighbor
        double inflow = rain_per_timestep;
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
                        terrain_elevations,
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
            terrain_elevations,
            elev_dimens,
            water_levels_in,
            grid_dimens,
            dt
        );
        double outflow = discharge.north + discharge.south + discharge.east + discharge.west;

        const double next_water =
            static_cast<double>(water_levels_in[pos.y * grid_dimens.x + pos.x]) + inflow - outflow;
        water_levels_out[pos.y * grid_dimens.x + pos.x] = static_cast<float>(next_water < 0.0 ? 0.0 : next_water);
    }
}
