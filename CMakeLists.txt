cmake_minimum_required(VERSION 3.21)
project(DLS)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#add_compile_options(-fprofile-generate)
#add_link_options(-fprofile-generate)

add_executable(DLS main.cpp trash.hpp DLS.hpp ComputationTask.hpp)
