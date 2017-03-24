//
// Created by lurker on 3/22/17.
//

#include "electric.h"

using namespace bbfmm;

void electric(Grid& g, levelset& ls, Surface& surf, Molecule& mol, scalar_t rescale){

#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());
#endif

    vector<point> source, target;
    vector<scalar_t > weight, normalX, normalY, normalZ;

    scalar_t dx = ls.dx / rescale;
    /*
     * map all points back to the actual protein surface.
     */
    for (index_t id = 0; id < surf.nodes.size(); ++id) {
        source.push_back(
                {
                    surf.nodes[id].data[0]/rescale,
                    surf.nodes[id].data[1]/rescale,
                    surf.nodes[id].data[2]/rescale
                }
        );

        target.push_back(
                {
                    surf.nodes[id].data[0]/rescale,
                    surf.nodes[id].data[1]/rescale,
                    surf.nodes[id].data[2]/rescale
                }
        );

        weight.push_back(surf.weight[id] * rescale * dx * SQR(dx));

        /*
         * normal vectors do not rescale.
         */
        normalX.push_back(surf.normals[id].data[0]);
        normalY.push_back(surf.normals[id].data[1]);
        normalZ.push_back(surf.normals[id].data[2]);

    }

    /*
     * setup for FMM.
     */
    index_t np = 8;
    index_t maxPoint = 320;
    index_t maxLevel = 10;

    scalar_t kappa = 0.;
    scalar_t dE = 80.0;
    scalar_t dI = 2.0;

    index_t  N = (index_t) source.size();

    scalar_t vacant_radius = 2 * dx;

    scalar_t area = std::accumulate(weight.begin(), weight.end(), 0.);
    std::cout << std::setw(15)<< "AREA" << " " << std::setw(8) << area << " A^2" <<std::endl;

    auto eval_G0 = [&](point& a, point& b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) return 0.;
        else {
            return 1.0 / 4.0 / M_PI / r;
        }

    };


    auto eval_Gk = [&](point& a, point& b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) return 0.;
        else {
            return exp(-kappa * r) / 4.0 / M_PI / r;
        }
    };

    auto eval_pG0x = [&](point& a, point& b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) return 0.;
        else return -(a.x  - b.x) / 4.0 / M_PI / r / d;
    };

    auto eval_pG0y = [&](point& a, point& b){
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) return 0.;
        else return -(a.y - b.y) / 4.0 / M_PI / r / d;
    };

    auto eval_pG0z = [&](point& a, point& b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) return 0.;
        else return -(a.z - b.z) / 4.0 / M_PI / r / d;
    };

    auto eval_pGkx = [&](point& a, point& b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) return 0.;
        else return -(a.x - b.x) * exp(-kappa * r) * (kappa * r + 1) / 4.0 / M_PI / r / d;
    };

    auto eval_pGky = [&](point& a, point& b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) return 0.;
        else return -(a.y - b.y) * exp(-kappa * r) * (kappa * r + 1) / 4.0 / M_PI / r / d;
    };

    auto eval_pGkz = [&](point& a, point& b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) return 0.;
        else return -(a.z - b.z) * exp(-kappa * r) * (kappa * r + 1) / 4.0 / M_PI / r / d;
    };

    auto maping = [&](vector<point>& _source, vector<point>& _target, vector<scalar_t>& _weight,
                      vector<scalar_t>& _normalX,  vector<scalar_t>& _normalY, vector<scalar_t>& _normalZ, Vector& _phi) {
        index_t _N = (index_t) _source.size();
        assert(_phi.row() == 2 * _N);
        kernel G0, Gk, pG0x, pG0y, pG0z, pGkx, pGky, pGkz;

        G0.eval = eval_G0; Gk.eval = eval_Gk;
        pG0x.eval = eval_pG0x; pG0y.eval = eval_pG0y;
        pG0z.eval = eval_pG0z; pGkx.eval = eval_pGkx;
        pGky.eval = eval_pGky; pGkz.eval = eval_pGkz;


        Vector pphi_pn(_N);
        Vector phiX(_N), phiY(_N), phiZ(_N);

        for (auto id = 0; id < _N; ++id) {
            pphi_pn(id) = _phi(id + _N) * _weight[id];
            phiX(id)    = _phi(id) * _normalX[id] * _weight[id];
            phiY(id)    = _phi(id) * _normalY[id] * _weight[id];
            phiZ(id)    = _phi(id) * _normalZ[id] * _weight[id];
        }

        G0.initialize(np, _source, _target, pphi_pn, _N, _N, maxPoint, maxLevel);
        pG0x.initialize(np, _source, _target, phiX, _N, _N, maxPoint, maxLevel);
        pG0y.initialize(np, _source, _target, phiY, _N, _N, maxPoint, maxLevel);
        pG0z.initialize(np, _source, _target, phiZ, _N, _N, maxPoint, maxLevel);

        Gk.initialize(np, _source, _target, pphi_pn, _N, _N, maxPoint, maxLevel);
        pGkx.initialize(np, _source, _target, phiX, _N, _N, maxPoint, maxLevel);
        pGky.initialize(np, _source, _target, phiY, _N, _N, maxPoint, maxLevel);
        pGkz.initialize(np, _source, _target, phiZ, _N, _N, maxPoint, maxLevel);

        Vector retG0, retGk, retG0X, retG0Y, retG0Z, retGkX, retGkY, retGkZ;
        G0.run(retG0); Gk.run(retGk);
        pG0x.run(retG0X); pG0y.run(retG0Y);
        pG0z.run(retG0Z); pGkx.run(retGkX);
        pGky.run(retGkY); pGkz.run(retGkZ);

        Vector output(2 * _N);

        for (auto id = 0; id < _N; ++id) {
            output(id) = 0.5 * _phi(id) + retG0X(id) + retG0Y(id) + retG0Z(id) - retG0(id);
            output(id + _N) = 0.5 * _phi(id) - retGkX(id) - retGkY(id) - retGkZ(id) + dI/dE * retGk(id);
        }

        return output;
    };

    auto FullMap = [&](Vector& phi) {
        return maping(source, target, weight, normalX, normalY, normalZ, phi);
    };


    Vector start(2 * N); setValue(start, 0.);

    Vector load(2 * N); setValue(load, 0.);
    for (auto id = 0; id < N; ++id) {
        for (auto atom_id = 0; atom_id <  mol.N; ++atom_id) {
            scalar_t d = SQR(source[id].x - mol.centers[atom_id].data[0]) +
                         SQR(source[id].y - mol.centers[atom_id].data[1]) +
                         SQR(source[id].z - mol.centers[atom_id].data[2]);

            scalar_t r = sqrt(d);
            load(id) += mol.charges[atom_id] / dI / 4.0 / M_PI / r;
        }
    }


    GMRES(FullMap, start, load, 20, 200, 1e-5);

    /*
     * solution is stored in start.
     */

}