/**
 * read-geo-data.cpp equires C++17 or greater.
 * Compile with `nvcc --std=c++17 cuda.cu -o cuda -lgdal -lstdc++fs`
 * after `module load cuda gdal`
 */

#include "kernels.cu"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include "read-geo-data.cpp"
#include "write-bmp.cpp"
#include "write-video.cpp"

int run_simulation(
    double *elevations_in,
    uint2 elev_dimens,
    double *water_levels_out,
    int num_timesteps,
    double dt,
    double total_rainfall_meters,
    struct WriteVideoContext *video_ctx
)
{
    auto log_cuda_error = [](const char *context, cudaError_t err) {
        fprintf(stderr, "%s; %s (ret=%d)\n", context, cudaGetErrorString(err), static_cast<int>(err));
    };

    int block_size_x = 32;
    int block_size_y = 32;

    uint2 grid_dimens = {
        .x = elev_dimens.x - 1,
        .y = elev_dimens.y - 1
    };

    // Add water each timestep to simulate rainfall
    double rain_per_timestep = total_rainfall_meters / num_timesteps;

    float *terrain_tile_d = nullptr;
    float *water_in_tile_d = nullptr;
    float *water_out_tile_d = nullptr;

    cudaError_t ret = cudaSuccess;
    int host_err = 0;
    size_t elev_count = static_cast<size_t>(elev_dimens.x) * static_cast<size_t>(elev_dimens.y);
    size_t grid_count = static_cast<size_t>(grid_dimens.x) * static_cast<size_t>(grid_dimens.y);

    std::vector<float> elevations_f(elev_count);
    dim3 threadsPerBlock(block_size_x, block_size_y);

    for (size_t i = 0; i < elev_count; i++) {
        elevations_f[i] = static_cast<float>(elevations_in[i]);
    }

    std::vector<float> water_levels_a_h(grid_count, 0.0f);
    std::vector<float> water_levels_b_h(grid_count, 0.0f);

    // Process in bounded tiles so GPU memory use is independent of full DEM size.
    const unsigned int tile_cells_x = 2048;
    const unsigned int tile_cells_y = 1024;
    const unsigned int max_tile_grid_x = tile_cells_x + 2;
    const unsigned int max_tile_grid_y = tile_cells_y + 2;
    const unsigned int max_tile_elev_x = max_tile_grid_x + 1;
    const unsigned int max_tile_elev_y = max_tile_grid_y + 1;
    const size_t max_tile_grid_count =
        static_cast<size_t>(max_tile_grid_x) * static_cast<size_t>(max_tile_grid_y);
    const size_t max_tile_elev_count =
        static_cast<size_t>(max_tile_elev_x) * static_cast<size_t>(max_tile_elev_y);

    std::vector<float> terrain_tile_h(max_tile_elev_count);
    std::vector<float> water_in_tile_h(max_tile_grid_count);
    std::vector<float> water_out_tile_h(max_tile_grid_count);
    const std::vector<float> *result_h = nullptr;

    ret = cudaMalloc((void **) &terrain_tile_d, max_tile_elev_count * sizeof(float));
    if (ret != cudaSuccess) {
        log_cuda_error("Device allocation failed for terrain tile buffer", ret);
        goto cleanup;
    }
    ret = cudaMalloc((void **) &water_in_tile_d, max_tile_grid_count * sizeof(float));
    if (ret != cudaSuccess) {
        log_cuda_error("Device allocation failed for water-in tile buffer", ret);
        goto cleanup;
    }
    ret = cudaMalloc((void **) &water_out_tile_d, max_tile_grid_count * sizeof(float));
    if (ret != cudaSuccess) {
        log_cuda_error("Device allocation failed for water-out tile buffer", ret);
        goto cleanup;
    }

    for (int i = 0; i < num_timesteps && host_err == 0; i++)
    {
        std::vector<float> &water_in_h = (i % 2 == 0) ? water_levels_a_h : water_levels_b_h;
        std::vector<float> &water_out_h = (i % 2 == 0) ? water_levels_b_h : water_levels_a_h;

        for (unsigned int y0 = 0; y0 < grid_dimens.y; y0 += tile_cells_y) {
            const unsigned int y1 = std::min(y0 + tile_cells_y, grid_dimens.y);
            const unsigned int interior_h = y1 - y0;
            const unsigned int halo_top = (y0 > 0) ? 1U : 0U;
            const unsigned int halo_bottom = (y1 < grid_dimens.y) ? 1U : 0U;

            for (unsigned int x0 = 0; x0 < grid_dimens.x; x0 += tile_cells_x) {
                const unsigned int x1 = std::min(x0 + tile_cells_x, grid_dimens.x);
                const unsigned int interior_w = x1 - x0;
                const unsigned int halo_left = (x0 > 0) ? 1U : 0U;
                const unsigned int halo_right = (x1 < grid_dimens.x) ? 1U : 0U;

                const unsigned int tile_origin_x = x0 - halo_left;
                const unsigned int tile_origin_y = y0 - halo_top;
                const unsigned int tile_grid_w = interior_w + halo_left + halo_right;
                const unsigned int tile_grid_h = interior_h + halo_top + halo_bottom;
                const unsigned int tile_elev_w = tile_grid_w + 1;
                const unsigned int tile_elev_h = tile_grid_h + 1;

                for (unsigned int row = 0; row < tile_elev_h; row++) {
                    const size_t src_idx =
                        static_cast<size_t>(tile_origin_y + row) * elev_dimens.x + tile_origin_x;
                    const size_t dst_idx = static_cast<size_t>(row) * tile_elev_w;
                    std::memcpy(
                        &terrain_tile_h[dst_idx],
                        &elevations_f[src_idx],
                        tile_elev_w * sizeof(float)
                    );
                }

                for (unsigned int row = 0; row < tile_grid_h; row++) {
                    const size_t src_idx =
                        static_cast<size_t>(tile_origin_y + row) * grid_dimens.x + tile_origin_x;
                    const size_t dst_idx = static_cast<size_t>(row) * tile_grid_w;
                    std::memcpy(
                        &water_in_tile_h[dst_idx],
                        &water_in_h[src_idx],
                        tile_grid_w * sizeof(float)
                    );
                }

                ret = cudaMemcpy(
                    terrain_tile_d,
                    terrain_tile_h.data(),
                    static_cast<size_t>(tile_elev_w) * tile_elev_h * sizeof(float),
                    cudaMemcpyHostToDevice
                );
                if (ret != cudaSuccess) {
                    log_cuda_error("Host to device memory copy failed for terrain tile", ret);
                    goto cleanup;
                }

                ret = cudaMemcpy(
                    water_in_tile_d,
                    water_in_tile_h.data(),
                    static_cast<size_t>(tile_grid_w) * tile_grid_h * sizeof(float),
                    cudaMemcpyHostToDevice
                );
                if (ret != cudaSuccess) {
                    log_cuda_error("Host to device memory copy failed for water tile", ret);
                    goto cleanup;
                }

                uint2 tile_elev_dimens = {tile_elev_w, tile_elev_h};
                uint2 tile_grid_dimens = {tile_grid_w, tile_grid_h};
                dim3 blocksPerCellGrid(
                    ceil(tile_grid_w / static_cast<double>(block_size_x)),
                    ceil(tile_grid_h / static_cast<double>(block_size_y))
                );

                timestep_forward_kernel<<<blocksPerCellGrid, threadsPerBlock>>>(
                    terrain_tile_d, tile_elev_dimens, water_in_tile_d, water_out_tile_d, tile_grid_dimens, dt, rain_per_timestep
                );

                ret = cudaDeviceSynchronize();
                if (ret != cudaSuccess) {
                    log_cuda_error("Kernel execution failed", ret);
                    goto cleanup;
                }

                ret = cudaMemcpy(
                    water_out_tile_h.data(),
                    water_out_tile_d,
                    static_cast<size_t>(tile_grid_w) * tile_grid_h * sizeof(float),
                    cudaMemcpyDeviceToHost
                );
                if (ret != cudaSuccess) {
                    log_cuda_error("Device to host memory copy failed for water tile", ret);
                    goto cleanup;
                }

                for (unsigned int row = 0; row < interior_h; row++) {
                    const size_t dst_idx = static_cast<size_t>(y0 + row) * grid_dimens.x + x0;
                    const size_t src_idx = static_cast<size_t>(row + halo_top) * tile_grid_w + halo_left;
                    std::memcpy(
                        &water_out_h[dst_idx],
                        &water_out_tile_h[src_idx],
                        interior_w * sizeof(float)
                    );
                }
            }
        }

        if (i % 10 == 0) {
            printf("Iteration %d\n", i);
        }

        if (video_ctx != nullptr) {
            if (write_video_frame(video_ctx, static_cast<int>(grid_dimens.x),
                                  static_cast<int>(grid_dimens.y),
                                  water_out_h.data()) != 0) {
                host_err = 1;
            }
        }
    }

    result_h = (num_timesteps % 2 == 0) ? &water_levels_a_h : &water_levels_b_h;

    for (size_t i = 0; i < grid_count; i++) {
        water_levels_out[i] = static_cast<double>((*result_h)[i]);
    }

cleanup:
    if (terrain_tile_d != nullptr) cudaFree(terrain_tile_d);
    if (water_in_tile_d != nullptr) cudaFree(water_in_tile_d);
    if (water_out_tile_d != nullptr) cudaFree(water_out_tile_d);

    if (ret != cudaSuccess)
        return EXIT_FAILURE;
    return host_err == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

void usage()
{
    std::cout << "USAGE: ./cuda <elev.tif> <num_iter> <total_rainfall(inches)>\n";
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        usage();
        return 1;
    }
    char *filename = argv[1];
    int num_timesteps = 1000; // defaults
    double rain_inches_total = 1;
    if (argc > 2)
    {
        num_timesteps = std::stoi(argv[2]);
    }
    if (argc > 3)
    {
        rain_inches_total = std::stod(argv[3]);
    }

    auto [elevations, elev_width, elev_height] = read_geo_data(filename);
    if (elevations == nullptr)
    {
        return EXIT_FAILURE;
    }
    unsigned int grid_width = (unsigned int) elev_width - 1;
    unsigned int grid_height = (unsigned int) elev_height - 1;
    uint2 elev_dimens = {(unsigned int) elev_width, (unsigned int) elev_height};

    printf("Loaded geo data; %dx%d\n", elev_width, elev_height);

    // Simulation Parameters
    double dt = 1; // Time step (seconds)

    double *water_levels = (double *) malloc(grid_width*grid_height*sizeof(double));

    const double fps = (dt > 0.0) ? (1.0 / dt) : 30.0;
    WriteVideoContext *video =
        write_video_begin(filename, static_cast<int>(grid_width), static_cast<int>(grid_height), fps);
    if (video == nullptr) {
        fprintf(stderr, "Failed to start ffmpeg video encoder (is ffmpeg on PATH?).\n");
        free(water_levels);
        return EXIT_FAILURE;
    }

    // Run the simulation
    int ret = run_simulation(
        elevations,
        elev_dimens,
        water_levels,
        num_timesteps,
        dt,
        inch_to_meter(rain_inches_total),
        video
    );
    if (write_video_end(video) != 0) {
        fprintf(stderr, "Video encoder finished with an error.\n");
        ret = EXIT_FAILURE;
    }
    if (ret != EXIT_SUCCESS) {
        free(water_levels);
        return EXIT_FAILURE;
    }

    // Cleanup
    free(water_levels);

    return 0;
}
