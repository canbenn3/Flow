#include <filesystem>
#include "gdal_priv.h"
#include <iostream>

namespace fs = std::filesystem;

struct GeoData
{
    double *elevations;
    int width;
    int height;
};

GeoData read_geo_data(const char *filename_read)
{
    GDALAllRegister();
    fs::path p(filename_read);
    GDALDataset *dataset = (GDALDataset *)GDALOpen(filename_read, GA_ReadOnly);

    GDALRasterBand *elevation_band = dataset->GetRasterBand(1);

    int width = elevation_band->GetXSize();
    int height = elevation_band->GetYSize();

    double *elevations = (double *)malloc(width * height * sizeof(double));

    CPLErr err = elevation_band->RasterIO(GF_Read, 0, 0, width, height, elevations, width, height, GDT_Float64, 0, 0);

    GDALClose(dataset);
    if (err != CE_None)
    {
        std::cerr << "Error reading data into the array.\n";
        return {NULL, 0, 0};
    }
    return {elevations, width, height};
}