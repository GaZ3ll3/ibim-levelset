# Implicit Boundary Integral Method with FMM

## setup files
If the system does not have ``glut.h``, then need to use CMakeLists_noview.txt.

- mv CMakeLists.txt CMakeLists.txt.bk
- mv CMakeLists_noview.txt CMakeLists.txt

## Compile
- mkdir build && cd build
- cmake .. 
- make

## Run it
There is an example config file in ``./data`` directory. Setting the numeric values and run with
- ./levelset ../data/input.cfg

