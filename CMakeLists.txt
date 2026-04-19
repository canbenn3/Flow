cmake_minimum_required(VERSION 3.18)

# Initialize project with both C++ and CUDA languages enabled
project(KinematicWaveSim LANGUAGES CXX CUDA)

# Require C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Required Packages
# find_package(OpenMP REQUIRED)
find_package(GDAL REQUIRED)

# ------------------------------------------------------------------------------
# 1. Serial Implementation
# ------------------------------------------------------------------------------
# add_executable(serial serial.cpp)
# target_include_directories(serial PRIVATE ${GDAL_INCLUDE_DIRS})
# Linking stdc++fs is often required for <filesystem> on older compilers (GCC < 9)
# target_link_libraries(serial PRIVATE ${GDAL_LIBRARIES} stdc++fs)

# ------------------------------------------------------------------------------
# 2. OpenMP Implementation
# ------------------------------------------------------------------------------
# add_executable(mp mp_impl.cpp)
# target_include_directories(mp PRIVATE ${GDAL_INCLUDE_DIRS})
# target_link_libraries(mp PRIVATE OpenMP::OpenMP_CXX ${GDAL_LIBRARIES} stdc++fs)

# ------------------------------------------------------------------------------
# 3. CUDA Implementation
# ------------------------------------------------------------------------------
add_executable(cuda cuda.cu)
# Use 'native' to compile for the GPU architecture currently installed on the build machine
set_target_properties(cuda PROPERTIES CUDA_ARCHITECTURES "native")
target_include_directories(cuda PRIVATE ${GDAL_INCLUDE_DIRS})
target_link_libraries(cuda PRIVATE ${GDAL_LIBRARIES} stdc++fs)