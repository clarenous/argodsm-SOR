cmake_minimum_required(VERSION 3.20)
project(sample_mpi)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

find_package(MPI)

include_directories(${MPI_INCLUDE_PATH})

add_executable(sample_mpi main.cpp)
