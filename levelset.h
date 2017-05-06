//
// Created by lurker on 2/23/17.
//

#ifndef LEVELSET_LEVELSET_H
#define LEVELSET_LEVELSET_H

#include "utils.h"
#include "molecule.h"
#include "Config.h"

class Grid {
public:
    index_t Nx;
    index_t Ny;
    index_t Nz;
    scalar_t * data;
    Grid(index_t _Nx, index_t _Ny, index_t _Nz);
    Grid(scalar_t val, index_t _Nx, index_t _Ny, index_t _Nz);
    ~Grid();
    scalar_t get(index_t i, index_t j, index_t k);
    void set(scalar_t val, index_t i, index_t j, index_t k);

    void output(std::string filename);
};

class levelset {
public:
    scalar_t dx;
    scalar_t sx;
    scalar_t sy;
    scalar_t sz;
    index_t Nx;
    index_t Ny;
    index_t Nz;
    index_t bandwidth;
    index_t shift;

    index_t thickness;
    scalar_t thres;

    levelset(index_t _Nx, index_t _Ny, index_t _Nz, index_t  _band, scalar_t _sx, scalar_t _sy,scalar_t _sz, scalar_t _dx, Config& cfg);

    ~levelset();
    void expand(Molecule& mol, Grid &g, scalar_t probe);
    void evolve(Grid& g, scalar_t final_t, scalar_t vel, scalar_t cfl_thres);
    void reinitialize(Grid &g, Grid &phi0, scalar_t final_t, scalar_t vel, scalar_t cfl_thres);

    void setWindow(Grid& g, scalar_t* window, index_t i, index_t j, index_t k);
    void setGradient(index_t dir, scalar_t* window, ls_point& uxp, ls_point& uxn);
    scalar_t getNorm(ls_point& Dun, ls_point& Dup);

    index_t countGradient(Grid& g, scalar_t thickness, scalar_t thres, scalar_t* window, bool display= false);
};

class Surface {
public:
    Surface(Grid& g, levelset& ls);

    ~Surface();
    vector<ls_point> nodes;
    vector<ls_point> normals;
    vector<scalar_t> weight;
};


#endif //LEVELSET_LEVELSET_H
