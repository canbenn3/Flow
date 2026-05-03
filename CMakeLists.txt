cmake_minimum_required(VERSION 3.10)
project(FlowSolver LANGUAGES CXX CUDA)

# 1. Set C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
find_package(GDAL REQUIRED)
find_package(MPI REQUIRED)

add_executable(mp 
    mp_impl.cpp 
    read-geo-data.cpp 
    write-bmp.cpp 
    write-jpg.cpp 
)
target_link_libraries(mp PRIVATE
    GDAL::GDAL
    OpenMP::OpenMP_CXX
    stdc++fs
)

add_executable(mpi
    mpi_impl.cpp 
    read-geo-data.cpp 
    write-bmp.cpp 
    write-jpg.cpp 
)
target_link_libraries(mpi PRIVATE
    GDAL::GDAL
    MPI::MPI_CXX
    stdc++fs
)

add_executable(serial 
    serial.cpp 
    read-geo-data.cpp 
    write-bmp.cpp 
    write-jpg.cpp 
)
target_link_libraries(serial PRIVATE GDAL::GDAL stdc++fs)

# 4. Single-GPU Cuda (Fixed naming from cuda -> cuda_sim)
add_executable(cuda_sim cuda.cu)
set_target_properties(cuda_sim PROPERTIES CUDA_ARCHITECTURES "native")
# Linking to GDAL::GDAL automatically handles include directories
target_link_libraries(cuda_sim PRIVATE GDAL::GDAL stdc++fs OpenMP::OpenMP_CXX)
target_compile_options(cuda_sim PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-O3 --use_fast_math --expt-relaxed-constexpr -Xcompiler=-fopenmp>
)

# 5. Distributed Cuda
# add_executable(cuda-distributed cuda-distributed.cu)

# Link MPI, GDAL, and the CUDA Runtime
# target_link_libraries(cuda-distributed 
#     PRIVATE 
#     MPI::MPI_CXX 
#     GDAL::GDAL
#     cudart
#     stdc++fs
# )

# target_compile_options(cuda-distributed PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
#     --generate-line-info
#     -O3
# >)

# Compiles support for Ampere (80), Ada Lovelace (89), and Hopper (90)
# set_target_properties(cuda-distributed PROPERTIES CUDA_ARCHITECTURES "80;89;90")
# target_include_directories(cuda-distributed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
