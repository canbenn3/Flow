/*
parse.cpp
FILE DESCRIPTION
Takes a given `.tif` representing GIS data and turns it into a topographical map with elevation being represented by a color scale
*/

#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include "gdal_priv.h"

#pragma pack(push, 1)
struct BMPFileHeader
{
    u_int16_t file_type{0x4D42};
    u_int32_t file_size{0};
    u_int16_t reserved1{0};
    u_int16_t reserved2{0};
    u_int32_t offset_data{54};
};
struct BMPInfoHeader
{
    u_int32_t size{40};
    int32_t width{0};
    int32_t height{0};
    u_int16_t planes{1};
    u_int16_t bit_count{24};
    u_int32_t compression{0};
    u_int32_t size_image{0};
    int32_t x_pixels_per_meter{0};
    int32_t y_pixels_per_meter{0};
    u_int32_t colors_used{0};
    u_int32_t colors_important{0};
};
#pragma pack(pop)

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    // 1. Register all GDAL format drivers
    GDALAllRegister();

    // 2. Open the GeoTIFF file
    const char *filename_r = argv[1];
    fs::path p(filename_r);
    std::string output_file = "output/" + p.stem().string() + ".bmp";
    std::cout << output_file << std::endl;
    GDALDataset *dataset = (GDALDataset *)GDALOpen(filename_r, GA_ReadOnly);

    // 3. Extract the Elevation Band
    // DEMs usually only have 1 band (unlike RGB photos which have 3)
    GDALRasterBand *elevationBand = dataset->GetRasterBand(1);

    // ... [Previous code opening the dataset remains the same] ...

    int fullWidth = elevationBand->GetXSize();
    int fullHeight = elevationBand->GetYSize();

    // 1. Define your N x N window size
    int nX = 1000; // Number of columns to read
    int nY = 1000; // Number of rows to read

    // Safety check: ensure N isn't larger than the actual map
    nX = std::min(nX, fullWidth);
    nY = std::min(nY, fullHeight);

    long long totalPixels = (long long)fullWidth * fullHeight;
    std::cout << "Reading a subset of " << nX << " x " << nY << " pixels.\n";

    // 2. Allocate memory for just the subset
    std::vector<float> elevationArray(totalPixels);

    // 3. Read only the N x N window
    CPLErr err = elevationBand->RasterIO(GF_Read,
                                         0, 0,                  // X and Y Offset (0,0 is top-left)
                                         fullWidth, fullHeight, // How many pixels to read from the file
                                         elevationArray.data(),
                                         fullWidth, fullHeight, // Size of the destination array in RAM
                                         GDT_Float32,
                                         0, 0);

    if (err != CE_None)
    {
        std::cerr << "Error reading data into the array." << std::endl;
        GDALClose(dataset);
        return 1;
    }

    std::cout << "Success! Data loaded into RAM.\n";

    // --- PROOF IT WORKS ---
    // In a 1D array, a 2D coordinate (x, y) is mapped as: index = (y * width) + x
    int testX = nX / 2;
    int testY = nY / 2;
    int index = (testY * nX) + testX;

    std::cout << "Elevation at the center of the map: "
              << elevationArray[index] << " meters.\n";

    std::cout << "converting elevation to color" << std::endl;

    const int width = fullWidth;
    const int height = fullHeight;

    BMPFileHeader fileHeader;
    BMPInfoHeader infoHeader;

    infoHeader.width = width;
    infoHeader.height = height;

    int paddingSize = (4 - (width * 3) % 4) % 4;
    u_int32_t fileSize = 54 + (3 * width + paddingSize) * height;
    fileHeader.file_size = fileSize;

    std::ofstream outFile(output_file, std::ios::binary);
    if (!outFile)
    {
        std::cerr << "Error: Could not open file for writing." << std::endl;
        return 1;
    }

    outFile.write(reinterpret_cast<const char *>(&fileHeader), sizeof(fileHeader));
    outFile.write(reinterpret_cast<const char *>(&infoHeader), sizeof(infoHeader));

    // TODO: getting min / max height takes a bit. See if can be parameterized in args
    float max_height = 0.0f;
    float min_height = 10000.0f;

    for (int i = 0; i < height * width; i++)
    {
        float factor = elevationArray[i];
        if (factor > max_height)
        {
            max_height = factor;
        }
        if (factor < min_height && factor > 0)
        {
            min_height = factor;
        }
    }

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float factor = static_cast<float>(elevationArray[(y * width) + x]);
            u_int16_t alt = static_cast<u_int16_t>(std::clamp(factor, min_height, max_height));
            u_int8_t high_byte = (alt >> 8) & 0xFF;
            u_int8_t low_byte = alt & 0xFF;
            u_int8_t blue_byte = 0;

            outFile.put(blue_byte);
            outFile.put(low_byte);
            outFile.put(high_byte);
        }
        for (int p = 0; p < paddingSize; ++p)
        {
            outFile.put(0);
        }
    }

    outFile.close();
    std::cout << "Success! Created " << output_file << " (" << width << "x" << height << ")" << std::endl;
    std::cout << "max height found: " << max_height << std::endl;
    std::cout << "min height found: " << min_height << std::endl;
    // 6. Clean up
    GDALClose(dataset);

    return 0;
}

void process_pixel(float &pixel)
{
    std::cout << "processing pixel" << std::endl;
}