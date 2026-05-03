/**
 * Encodes a simulation visualization as H.264/MP4 by piping raw RGB24 frames to ffmpeg.
 * Pixel mapping matches write-jpg.cpp (blue channel = water depth in inches, clamped).
 * Requires `ffmpeg` on PATH at runtime; RGB packing uses OpenMP on the CPU.
 */

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include <sys/wait.h>

double meter_to_inch(double meter);

namespace fs = std::filesystem;

struct WriteVideoContext
{
    FILE *pipe = nullptr;
    std::vector<uint8_t> scratch;
    /// Simulation grid size (water buffer is grid_w * grid_h).
    int grid_w = 0;
    int grid_h = 0;
    /// Padded to even dimensions for H.264 / yuv420p (>= grid, each side even).
    int enc_w = 0;
    int enc_h = 0;
};

namespace {

constexpr double kMaxClampInches = 12.0;

std::string shell_single_quote(const std::string &s)
{
    std::string out = "'";
    for (char c : s)
    {
        if (c == '\'')
            out += "'\\''";
        else
            out += c;
    }
    out += '\'';
    return out;
}

template <typename T>
int write_video_frame_impl(
    WriteVideoContext *ctx,
    int grid_width,
    int grid_height,
    const T *water_levels
)
{
    if (ctx == nullptr || ctx->pipe == nullptr)
        return 1;
    if (grid_width != ctx->grid_w || grid_height != ctx->grid_h)
        return 1;

    const int ew = ctx->enc_w;
    const int eh = ctx->enc_h;
    const int gw = ctx->grid_w;
    const int gh = ctx->grid_h;

    const size_t plane = static_cast<size_t>(ew) * static_cast<size_t>(eh) * 3u;
    if (ctx->scratch.size() != plane)
        return 1;

    uint8_t *dst = ctx->scratch.data();

#pragma omp parallel for schedule(static)
    for (int y = 0; y < eh; ++y)
    {
        for (int x = 0; x < ew; ++x)
        {
            const size_t base =
                (static_cast<size_t>(y) * static_cast<size_t>(ew) + static_cast<size_t>(x)) * 3u;
            if (x < gw && y < gh)
            {
                const size_t row = static_cast<size_t>(y) * static_cast<size_t>(gw);
                double water_level = meter_to_inch(
                    static_cast<double>(water_levels[row + static_cast<size_t>(x)]));
                water_level = std::max(0.0, std::min(water_level, kMaxClampInches));
                uint8_t intensity =
                    static_cast<uint8_t>((water_level / kMaxClampInches) * 255.0);
                dst[base + 0] = 0;
                dst[base + 1] = 0;
                dst[base + 2] = intensity;
            }
            else
            {
                dst[base + 0] = 0;
                dst[base + 1] = 0;
                dst[base + 2] = 0;
            }
        }
    }

    size_t n = std::fwrite(dst, 1, plane, ctx->pipe);
    return n == plane ? 0 : 1;
}

} // namespace

WriteVideoContext *write_video_begin(const char *tif_filename, int width, int height, double fps)
{
    if (width <= 0 || height <= 0 || fps <= 0.0)
        return nullptr;

    // libx264 + yuv420p needs even width and height; pad with black if needed.
    const int enc_w = (width + 1) & ~1;
    const int enc_h = (height + 1) & ~1;

    fs::create_directories("output");

    fs::path p(tif_filename);
    std::string output_file = "output/" + p.stem().string() + ".mp4";

    std::string cmd =
        "ffmpeg -hide_banner -loglevel error -y "
        "-f rawvideo -pixel_format rgb24 "
        "-video_size " + std::to_string(enc_w) + "x" + std::to_string(enc_h) + " "
        "-framerate " + std::to_string(fps) + " "
        "-i - "
        "-c:v libx264 -preset medium -crf 23 -pix_fmt yuv420p "
        + shell_single_quote(output_file);

    FILE *pipe = popen(cmd.c_str(), "w");
    if (pipe == nullptr)
        return nullptr;

    auto *ctx = new WriteVideoContext();
    ctx->pipe = pipe;
    ctx->grid_w = width;
    ctx->grid_h = height;
    ctx->enc_w = enc_w;
    ctx->enc_h = enc_h;
    ctx->scratch.resize(static_cast<size_t>(enc_w) * static_cast<size_t>(enc_h) * 3u);
    return ctx;
}

int write_video_frame(WriteVideoContext *ctx, int width, int height, double *water_levels)
{
    return write_video_frame_impl(ctx, width, height, water_levels);
}

int write_video_frame(WriteVideoContext *ctx, int width, int height, const float *water_levels)
{
    return write_video_frame_impl(ctx, width, height, water_levels);
}

int write_video_end(WriteVideoContext *ctx)
{
    if (ctx == nullptr)
        return 1;

    int status = 0;
    if (ctx->pipe != nullptr)
    {
        status = pclose(ctx->pipe);
        ctx->pipe = nullptr;
    }
    ctx->scratch.clear();
    ctx->scratch.shrink_to_fit();
    delete ctx;

    if (status == -1)
        return 1;
    if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
        return 0;
    return 1;
}
