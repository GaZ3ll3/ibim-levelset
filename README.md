# Implicit Boundary Integral Method with FMM

## OpenBlas
- OpenBlas has to be compiled with ``OPENMP=1``.
## Compile
- mkdir build && cd build
- cmake .. 
- make

## Run it
There is an example config file in ``./data`` directory. Setting the numeric values and run with
- ./levelset ../data/input.cfg

