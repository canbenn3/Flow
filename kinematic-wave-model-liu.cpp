/**
 * Implementation of the surface flow model proposed in:
 * 
 * Liu, Q., Chen, L., Li, J., & Singh, V. (2003). "Two-dimensional kinematic
 * wave model of overland-flow." Journal of Hydrology. 291(28–41).
 * https://doi.org/10.1016/j.jhydrol.2003.12.023
 */

#include <cmath>

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

GridCell flow_vector_direction(int northwest, int northeast, int southeast, int southwest) {
    // See Figure 3b for corner mappings
    int z1 = northwest;
    int z2 = northeast;
    int z3 = southeast;
    int z4 = southwest;

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

Discharge compute_discharge(GridCell cell, double cell_water_level) {
    // Very naive implementation simply scales outflow based on cell grade to
    // compute total cell discharge (q)
    // TODO: consider adding infiltration to the model
    double q = cell.theta/(pi/2.0)*cell_water_level;

    double qx = q * cos(cell.gamma);
    double qy = q * sin(cell.gamma);

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
