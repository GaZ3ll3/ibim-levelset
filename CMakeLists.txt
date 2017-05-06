cmake_minimum_required(VERSION 2.8)
project(levelset)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -DRUN_OMP -march=native -std=c++11 -fopenmp ")

find_package(GLUT)

if (GLUT_FOUND)
    set(SOURCE_FILES main.cpp levelset.cpp levelset.h utils.h molecule.h ls_point.h view.h view.cpp linalg.h blas_wraper.h gmres.h bbfmm.h bbfmm.cpp blas_wraper.cpp electric.cpp gmres.cpp Config.cpp Config.h)
    add_executable(levelset ${SOURCE_FILES})
    target_link_libraries(levelset m blas GL GLU glut)
else ()
    set(SOURCE_FILES main.cpp levelset.cpp levelset.h utils.h molecule.h ls_point.h linalg.h blas_wraper.h gmres.h bbfmm.h bbfmm.cpp blas_wraper.cpp electric.cpp gmres.cpp Config.cpp Config.h)
    add_executable(levelset ${SOURCE_FILES})
    target_link_libraries(levelset m blas)
endif ()

