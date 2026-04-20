cmake_minimum_required(VERSION 3.10)
project(FlowSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)
find_package(GDAL REQUIRED)
find_package(MPI REQUIRED)

add_executable(mp 
    mp_impl.cpp 
    read-geo-data.cpp 
    write-bmp.cpp 
    write-jpg.cpp 
    util.cpp
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
    util.cpp
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
    util.cpp
)
target_link_libraries(serial PRIVATE
    GDAL::GDAL
    stdc++fs
)