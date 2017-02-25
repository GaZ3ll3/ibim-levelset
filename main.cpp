#define DISP
//#define VIEW
#include <iostream>
#include "levelset.h"
#include "view.h"

int main() {

    Molecule mol; mol.load("../data/1hje.pqr");

    scalar_t s = mol.centralize(200.0); mol.getCenter();
    scalar_t pr = 1.4 * s;

    scalar_t grid_lo = -300.0, grid_hi = 300;

    index_t size = 256;
    scalar_t dx = (grid_hi - grid_lo) / scalar_t(size);

    levelset ls(size, size, size, 8, grid_lo, grid_lo, grid_lo, dx);

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

    RUN("REINIT", ls.reinitialize(g, phi0, 400, 1, 0.5));




#ifdef VIEW
    view& v = view::getInstance(2);
    v.loadLevelSet(ls, g);
    v.run();
#endif

    return 0;
}