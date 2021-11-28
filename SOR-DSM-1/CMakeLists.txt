cmake_minimum_required(VERSION 3.20)
project(sor_argo_1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_C_COMPILER mpicc)

include_directories(${ARGO_INSTALL_DIRECTORY}/include)
link_directories(${ARGO_INSTALL_DIRECTORY}/lib)

add_executable(sor_argo_1 sor_main.cc ../src/dsm-C.cc)
target_link_libraries(sor_argo_1 argo argobackend-mpi)
