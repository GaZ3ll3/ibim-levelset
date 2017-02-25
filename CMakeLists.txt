cmake_minimum_required(VERSION 2.8)
project(levelset)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -DRUN_OMP -march=native -std=c++11 -fopenmp ")

set(SOURCE_FILES main.cpp levelset.cpp levelset.h utils.h molecule.h point.h view.h view.cpp)
add_executable(levelset ${SOURCE_FILES})
target_link_libraries(levelset m blas GL GLU glut)