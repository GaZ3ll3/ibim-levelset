#define DISP
#define VIEW
#define GRID
#include <iostream>
#include "electric.h"

#undef VIEW
#undef GRID

#ifdef VIEW
#include "view.h"
#endif

int main(int argc, char* argv[]) {

#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());
#endif

    if (argc <= 1) {
        std::cout << "USE " <<argv[0] << " PATH_OF_CONFIG_FILE " << std::endl;
        exit(0);
    }

    Config cfg;
    std::ifstream cfgFile;cfgFile.open(argv[1], std::ifstream::in);
    cfg.parse(cfgFile);cfgFile.close();

    cfg.print();

    Molecule mol; mol.load(cfg.options["pqr_file"]);

    scalar_t s = mol.centralize(240.0);
    mol.getCenter();
    scalar_t pr = 1.4 * s;

    std::cout << std::setw(15) << "1 Angstrom" << " " << std::setw(8) << s << " grids" << std::endl;

    scalar_t grid_lo = -300.0, grid_hi = 300;
    index_t size = atoi(cfg.options["grid_size"].c_str());
    scalar_t dx = (grid_hi - grid_lo) / scalar_t(size);

    levelset ls(size, size, size, 8, grid_lo, grid_lo, grid_lo, dx, cfg);

    Grid g(-2.0 * size, ls.Nx, ls.Ny, ls.Nz);
    Grid phi0(0., ls.Nx, ls.Ny, ls.Nz);

    RUN("EXPAND", ls.expand(mol, g, pr));
    RUN("INWARD", ls.evolve(g, 1.4, s, 0.2));

    g *= -1.0;
    phi0 = g;

    RUN("REINIT 1st", ls.reinitialize(g, phi0, atoi(cfg.options["reinit_step"].c_str()), 1, 0.5));

    phi0 = g;

    RUN("REINIT 2nd", ls.reinitialize(g, phi0, atoi(cfg.options["reinit_step"].c_str()), 1, 0.5));

    Surface surf(g, ls, s);


#ifdef GRID
    g.output("../data/output.grid");
    surf.output("../data/output.node");
#endif


#ifdef VIEW
    view& v = view::getInstance(2);
    v.loadLevelSet(g, ls, surf);
    v.run();
#endif

    electric(g, ls, surf, mol, s, cfg);

    return 0;
}