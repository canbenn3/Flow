/*
generate_test.cpp
Used to generate testing dummy data for getting started.
*/
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include "gdal_priv.h"

const int SIZE = 5;

// Helper to write a 1D vector to a GeoTIFF
bool write_geotiff(const std::string &filename, const std::vector<float> &data)
{
    GDALDriver *driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!driver)
        return false;

    GDALDataset *dataset = driver->Create(filename.c_str(), SIZE, SIZE, 1, GDT_Float32, NULL);
    if (!dataset)
        return false;

    double geoTransform[6] = {0.0, 1.0, 0.0, 0.0, 0.0, -1.0};
    dataset->SetGeoTransform(geoTransform);

    GDALRasterBand *band = dataset->GetRasterBand(1);
    CPLErr err = band->RasterIO(GF_Write, 0, 0, SIZE, SIZE,
                                (void *)data.data(), SIZE, SIZE,
                                GDT_Float32, 0, 0);

    GDALClose(dataset);
    return err == CE_None;
}

int main()
{
    GDALAllRegister();

    // 1. Hill (Highest at center, sloping out)
    std::vector<float> hill(SIZE * SIZE);
    for (int i = 0; i < 25; i++)
    {
        int r = i / SIZE, c = i % SIZE;
        hill[i] = 20.0f - (std::abs(2 - r) + std::abs(2 - c));
    }

    // 2. Ridge (High horizontal line in middle, slopes North/South)
    std::vector<float> ridge(SIZE * SIZE);
    for (int i = 0; i < 25; i++)
    {
        int r = i / SIZE;
        ridge[i] = 15.0f - std::abs(2 - r);
    }

    // 3. Crest (Sharp peak at the very top row, sloping all the way down)
    std::vector<float> crest(SIZE * SIZE);
    for (int i = 0; i < 25; i++)
    {
        int r = i / SIZE;
        crest[i] = 10.0f - r;
    }

    // 4. Diagonal Crest (High at Top-Left [0,0], lowest at Bottom-Right [4,4])
    std::vector<float> diag_crest(SIZE * SIZE);
    for (int i = 0; i < 25; i++)
    {
        int r = i / SIZE, c = i % SIZE;
        diag_crest[i] = 10.0f - (r + c);
    }

    // 5. Diagonal Hill (A diagonal "spine" from Top-Left to Bottom-Right)
    std::vector<float> diag_hill(SIZE * SIZE);
    for (int i = 0; i < 25; i++)
    {
        int r = i / SIZE, c = i % SIZE;
        // Distance from the line r = c
        diag_hill[i] = 10.0f - std::abs(r - c);
    }

    std::vector<float> bowl(SIZE * SIZE);
    for (int i = 0; i < 25; i++)
    {
        int r = i / SIZE, c = i % SIZE;
        bowl[i] = (std::abs(2 - r) + std::abs(2 - c));
    }

    std::vector<float> valley(SIZE * SIZE);
    for (int i = 0; i < 25; i++)
    {
        int c = i % SIZE;
        valley[i] = std::abs(2 - c);
    }

    // Export all test cases
    struct TestCase
    {
        std::string name;
        std::vector<float> data;
    };
    std::vector<TestCase> tests = {
        {"hill.tif", hill},
        {"ridge.tif", ridge},
        {"crest.tif", crest},
        {"diag_crest.tif", diag_crest},
        {"diag_hill.tif", diag_hill},
        {"bowl.tif", bowl},
        {"valley.tif", valley}};

    for (const auto &t : tests)
    {
        std::string path = "./test_data/" + t.name;
        if (write_geotiff(path, t.data))
        {
            std::cout << "Generated: " << path << std::endl;
        }
        else
        {
            std::cerr << "Failed to generate: " << path << std::endl;
        }
    }

    return 0;
}