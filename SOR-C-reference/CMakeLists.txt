cmake_minimum_required(VERSION 3.20)
project(sor_c)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_COMPILER gcc)
set(CMAKE_C_COMPILER g++)

add_executable(sor_c sor_main.cc)
target_link_libraries(sor_c stdc++)
