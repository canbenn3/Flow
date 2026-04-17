#include <filesystem>
#include <iostream>
#include <fstream>
#include <cstdint>

namespace fs = std::filesystem;

double meter_to_inch(double meter);

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

// Writes a bitmap file (.bmp) into an output directory
int write_bmp(const char *tif_filename, int width, int height, double *water_levels)
{
    fs::path p(tif_filename);
    std::string output_file = "output/" + p.stem().string() + ".bmp";
    std::cout << "Writing bmp file to " << output_file << "\n";

    BMPFileHeader file_header;
    BMPInfoHeader info_header;
    info_header.width = width;
    info_header.height = height;

    // BMP rows must be a multiple of 4 bytes.
    int padding_size = (4 - (width * 3) % 4) % 4;
    uint32_t file_size = 54 + (3 * width + padding_size) * height;
    file_header.file_size = file_size;

    std::ofstream out_file(output_file, std::ios::binary);
    if (!out_file)
    {
        std::cerr << "Error: Could not open file for writing." << "\n";
        return 1;
    }

    out_file.write(reinterpret_cast<const char *>(&file_header), sizeof(file_header));
    out_file.write(reinterpret_cast<const char *>(&info_header), sizeof(info_header));

    // Save water levels to bmp. Consider 12 inches as the deepest blue and 0 inches is blank (black)
    double max_clamp = 12.0;

    for (int y = 0; y < height; ++y)
    {
        // BMP files are written bottom-to-top. Read from the array bottom-up
        // to prevent the final image from being flipped vertically.
        int array_y = (height - 1) - y;

        for (int x = 0; x < width; ++x)
        {
            double water_level = meter_to_inch(water_levels[(array_y * width) + x]);

            // Clamp the water level between 0 and 12 inches
            if (water_level < 0.0)
                water_level = 0.0;
            if (water_level > max_clamp)
                water_level = max_clamp;

            // Scale to a 0-255 byte value
            uint8_t blue_intensity = static_cast<uint8_t>((water_level / max_clamp) * 255.0);

            // BMP pixels are BGR (Blue, Green, Red)
            uint8_t pixel[3] = {blue_intensity, 0, 0};

            out_file.write(reinterpret_cast<const char *>(pixel), 3);
        }

        // Write the padding at the end of each row
        uint8_t padding[3] = {0, 0, 0};
        if (padding_size > 0)
        {
            out_file.write(reinterpret_cast<const char *>(padding), padding_size);
        }
    }

    out_file.close();
    std::cout << "Successfully wrote output file.\n";
    return 0;
}