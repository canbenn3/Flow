#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <filesystem>
#include <vector>
// #include "util.cpp"

double meter_to_inch(double meter);

namespace fs = std::filesystem;

int write_jpg(const char *tif_filename, int width, int height, double *water_levels)
{
    fs::path p(tif_filename);
    std::string output_file = "output/" + p.stem().string() + ".jpg";

    // 1. Create a buffer of 8-bit RGB pixels (stb doesn't need padding)
    std::vector<uint8_t> image_data(width * height * 3);
    double max_clamp = 12.0;

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            double water_level = meter_to_inch(water_levels[y * width + x]);

            // Clamp and Scale
            water_level = std::max(0.0, std::min(water_level, max_clamp));
            uint8_t intensity = static_cast<uint8_t>((water_level / max_clamp) * 255.0);

            int index = (y * width + x) * 3;
            image_data[index + 0] = 0;         // Red
            image_data[index + 1] = 0;         // Green
            image_data[index + 2] = intensity; // Blue
        }
    }

    // 2. Write the JPG (Quality: 90)
    int success = stbi_write_jpg(output_file.c_str(), width, height, 3, image_data.data(), 90);

    return success ? 0 : 1;
}