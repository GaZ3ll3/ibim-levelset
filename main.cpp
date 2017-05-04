#define DISP
//#define VIEW
#include <iostream>
#include "electric.h"

#ifdef VIEW
#include "view.h"
#endif

int main() {

    Config cfg;
    std::ifstream cfgFile;cfgFile.open("../data/input.cfg", std::ifstream::in);
    cfg.parse(cfgFile);cfgFile.close();

    Molecule mol; mol.load(cfg.options["pqr_file"]);

    scalar_t s = mol.centralize(200.0); mol.getCenter();
    scalar_t pr = 1.4 * s;
    scalar_t grid_lo = -300.0, grid_hi = 300;
    index_t size = atoi(cfg.options["grid_size"].c_str());
    scalar_t dx = (grid_hi - grid_lo) / scalar_t(size);

    levelset ls(size, size, size, 8, grid_lo, grid_lo, grid_lo, dx, cfg);

    Grid g(-2.0 * size, ls.Nx, ls.Ny, ls.Nz);
    Grid phi0(0., ls.Nx, ls.Ny, ls.Nz);

    RUN("EXPAND", ls.expand(mol, g, pr));
    RUN("INWARD", ls.evolve(g, 1.4, s, 0.5));

    for (index_t i = 0; i < ls.Nx; ++i) {
        for (index_t j = 0; j < ls.Ny; ++j) {
            for (index_t k = 0; k < ls.Nz; ++k) {
                index_t  I = i * ls.Ny * ls.Nz + j * ls.Nz + k;
                g.data[I] = -g.data[I];
                phi0.data[I] = g.data[I];
            }
        }
    }



    RUN("REINIT", ls.reinitialize(g, phi0, atoi(cfg.options["reinit_step"].c_str()), 1, 0.8));


    for (index_t i = 0; i < ls.Nx; ++i) {
        for (index_t j = 0; j < ls.Ny; ++j) {
            for (index_t k = 0; k < ls.Nz; ++k) {
                index_t  I = i * ls.Ny * ls.Nz + j * ls.Nz + k;
                phi0.data[I] = g.data[I];
            }
        }
    }

    RUN("REINIT", ls.reinitialize(g, phi0, atoi(cfg.options["reinit_step"].c_str()), 1, 0.8));

    /*
     * another run for inclusion removal.
     *
     * only to locate exterior point by marching from corner. And then
     * march the whole tube to take the tube out.
     *
     * use BFS.
     */
    

    Surface surf(g, ls);

//    std::ofstream gridFile;
//    gridFile.open("../data/test.grid"+ std::to_string(g.Nx) +std::to_string(g.Ny)+std::to_string(g.Nz));
//
//    for (index_t i = 0; i < ls.Nx; ++i) {
//        for (index_t j = 0; j < ls.Ny; ++j) {
//            for (index_t k = 0; k < ls.Nz; ++k) {
//                index_t  I = i * ls.Ny * ls.Nz + j * ls.Nz + k;
//                gridFile << g.data[I] << "\n";
//            }
//        }
//    }
//    gridFile.close();


#ifdef VIEW
    view& v = view::getInstance(2);
    v.loadLevelSet(g, ls, surf);
    v.run();
#endif

    electric(g, ls, surf, mol, s, cfg);




    return 0;
}