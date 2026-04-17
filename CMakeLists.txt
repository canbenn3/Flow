cmake_minimum_required(VERSION 3.10)
project(FlowSolver LANGUAGES CXX)

# 1. Set C++17 standard (Fixes the 'filesystem' and 'structured bindings' errors)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 2. Find Dependencies
find_package(OpenMP REQUIRED)
find_package(GDAL REQUIRED)

# 3. Define the executable
# List all your .cpp files here. This avoids the "redefinition" errors 
# by letting CMake handle the object linking properly.
add_executable(mp 
    mp_impl.cpp 
    read-geo-data.cpp 
    write-bmp.cpp 
    write-jpg.cpp 
    util.cpp
)

# 4. Link Libraries and Include Directories
target_link_libraries(mp PRIVATE GDAL::GDAL OpenMP::OpenMP_CXX stdc++fs)

# (Optional) Repeat for the serial version
add_executable(serial 
    serial.cpp 
    read-geo-data.cpp 
    write-bmp.cpp 
    write-jpg.cpp 
    util.cpp
)
target_link_libraries(serial PRIVATE GDAL::GDAL stdc++fs)