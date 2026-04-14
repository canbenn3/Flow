/**
 * Implementation of the surface flow model proposed in:
 * 
 * Liu, Q., Chen, L., Li, J., & Singh, V. (2003). "Two-dimensional kinematic
 * wave model of overland-flow." Journal of Hydrology. 291(28–41).
 * https://doi.org/10.1016/j.jhydrol.2003.12.023
 */

#include <cmath>
#include <algorithm> 

const double pi = 3.14159265358979323846;

struct GridCell {
    /// Grade (steepness) in radians
    double theta;
    /// Direction of flow (0=South, pi/2=East, etc.)
    double gamma;
};

struct Discharge {
    double north = 0.0;
    double east = 0.0;
    double south = 0.0;
    double west = 0.0;
};

GridCell flow_vector_direction(double northwest, double northeast, double southeast, double southwest) {
    // See Figure 3b for corner mappings
    double z1 = northwest;
    double z2 = northeast;
    double z3 = southeast;
    double z4 = southwest;

    // Equation 13
    double dF_dx_transformed = (z1+z2-z3-z4)/4.0;

    // Equation 14
    double dF_dy_transformed = (z1-z2-z3+z4)/4.0;

    // For the following equations:
    // Grid dimensions (dx and dy; defined in Figure 3a) are taken to be 1:1

    // Equation 11
    double theta = atan(sqrt(pow(dF_dx_transformed, 2)*4.0 + pow(dF_dy_transformed, 2)*4.0));

    // Equation 12
    double gamma = atan2(dF_dy_transformed, dF_dx_transformed);

    return {.theta = theta, .gamma = gamma};
}

// 'dt' (time step) scales flow rates into volumes
Discharge compute_discharge(GridCell cell, double cell_water_level, double dt) {
    // Calculate flow rate, then multiply by the time step (dt) to get the volume 
    // of water attempting to move during this iteration.
    double q_rate = (cell.theta / (pi / 2.0)) * cell_water_level;
    double q_vol = q_rate * dt;

    double qx = q_vol * cos(cell.gamma);
    double qy = q_vol * sin(cell.gamma);

    // --- MASS CONSERVATION CHECK ---
    // The total water trying to leave the cell is the sum of the absolute directional flows.
    double total_outflow = std::abs(qx) + std::abs(qy);

    // If the gradient is so steep that it tries to push out more water than exists,
    // we scale the vectors down proportionally. This prevents "negative water" bugs.
    if (total_outflow > cell_water_level) {
        double scale = cell_water_level / total_outflow;
        qx *= scale;
        qy *= scale;
    }

    Discharge discharge = {};

    if (qx > 0) {
        discharge.south = qx;
    } else {
        discharge.north = -qx;
    }

    if (qy > 0) {
        discharge.east = qy;
    } else {
        discharge.west = -qy;
    }

    return discharge;
}