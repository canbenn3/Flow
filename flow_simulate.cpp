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
struct GridCell
{
    float elevation;
    float surface_water;
    float absorbed_water;
    float absorption_rate;
    float total_elev() {
        return elevation + surface_water;
    };
    void do_absorb() {
        
    }
};
#pragma pack(pop)

namespace fs = std::filesystem;

void usage()
{
    std::cout << R"(
./flow_sim <file/path/*.tif> <rain_inches> <time_span(minutes)> [max_height] [min_height]
    )" << std::endl;
}

int main(int argc, char *argv[])
{
    GDALAllRegister();

    std::cout << "argc: " << argc << std::endl;
    if (argc < 4)
    {
        usage();
        return 1;
    }
    const char *filename_r = argv[1];
    const float rain_inches = std::stof(argv[2]);
    const float time_span = std::stof(argv[3]);
    float min_h = 1.0;
    float max_h = 10.0;
    if (argc > 4)
    {
        max_h = std::stof(argv[4]);
    }
    if (argc > 5)
    {
        min_h = std::stof(argv[5]);
    }

    std::cout << "doing flow sim of " << filename_r << " with min height " << min_h << " and max height " << max_h << std::endl;
    fs::path p(filename_r);
    std::string output_file = "output/" + p.stem().string() + ".bmp";
    GDALDataset *dataset = (GDALDataset *)GDALOpen(filename_r, GA_ReadOnly);

    if (dataset == nullptr)
    {
        std::cerr << "Error: Could not open the DEM file." << std::endl;
        return 1;
    }

    GDALRasterBand *elevationBand = dataset->GetRasterBand(1);
    int fullWidth = elevationBand->GetXSize();
    int fullHeight = elevationBand->GetYSize();
    long long totalPixels = (long long)fullWidth * fullHeight;

    std::vector<float> elevationArray(totalPixels);

    CPLErr err = elevationBand->RasterIO(GF_Read,
                                         0, 0,
                                         fullWidth, fullHeight,
                                         elevationArray.data(), fullWidth,
                                         fullHeight, GDT_Float32,
                                         0, 0);

    if (err != CE_None)
    {
        std::cerr << "Error reading data into the array." << std::endl;
        GDALClose(dataset);
        return 1;
    }

    BMPFileHeader fileHeader;
    BMPInfoHeader infoHeader;

    infoHeader.width = fullWidth;
    infoHeader.height = fullHeight;
    int paddingSize = (4 - (fullWidth * 3) % 4) % 4;
    u_int32_t fileSize = 54 + (3 * fullWidth + paddingSize) * fullHeight;
    std::ofstream outFile(output_file, std::ios::binary);
    if (!outFile)
    {
        std::cerr << "Error: Could not open file for writing." << std::endl;
        return 1;
    }

    outFile.write(reinterpret_cast<const char *>(&fileHeader), sizeof(fileHeader));
    outFile.write(reinterpret_cast<const char *>(&infoHeader), sizeof(infoHeader));
}