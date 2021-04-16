cmake_minimum_required(VERSION 3.1.0)
project(image_convolution)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -O3")
find_package(OpenMP)

add_executable(convolution pbm_transport.cpp)
target_link_libraries(convolution OpenMP::OpenMP_CXX)
