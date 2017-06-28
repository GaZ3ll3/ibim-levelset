//
// Created by lurker on 3/22/17.
//

#include "electric.h"

using namespace bbfmm;

void electric(Grid& g, levelset& ls, Surface& surf, Molecule& mol, scalar_t rescale, Config& cfg){
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
    std::cout << std::setw(15) <<"POINTS NUM"  << " " << std::setw(8) << source.size() << std::endl;

    /*
     * setup for FMM.
     */
    index_t np = atoi(cfg.options["fmm_np"].c_str());
    index_t maxPoint = atoi(cfg.options["fmm_maxpoint"].c_str());
    index_t maxLevel = atoi(cfg.options["fmm_maxlevel"].c_str());

    scalar_t kappa =atof(cfg.options["solvent_kappa"].c_str());
    scalar_t dE = atof(cfg.options["solvent_dE"].c_str());
    scalar_t dI = atof(cfg.options["solvent_dI"].c_str());

    index_t  N = (index_t) source.size();

    scalar_t vacant_radius = atof(cfg.options["tau"].c_str()) * dx;

    scalar_t area = std::accumulate(weight.begin(), weight.end(), 0.);
    std::cout << std::setw(15)<< "AREA APPROX" << " " << std::setw(8)<< area << " A^2" <<std::fixed<<std::endl;



    auto eval_G0 = [&](point& a, point& b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) {
            return 1.0 / 2.0 / M_PI / vacant_radius;
        }
        else {
            return 1.0 / 4.0 / M_PI / r;
        }

    };


    auto eval_Gk = [&](point& a, point& b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) {
            if (kappa == 0.) return 1.0/4.0/M_PI/vacant_radius;
            return (1 - exp(-kappa * vacant_radius)) / kappa / 2.0 / M_PI / vacant_radius / vacant_radius;
        }
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

/*
 *  dBIE, singular part can be regularized as anything, e.g. 0. First order is observed.
 */
    auto K1x = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        if (r < vacant_radius) return 0.;
        else return (a.x - b.x) * (dE / dI * exp(-kappa * r) * (kappa * r + 1) - 1.) / 4.0 / M_PI / r / d;
    };

    auto K1y = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        if (r < vacant_radius) return 0.;
        else return (a.y - b.y) * (dE / dI * exp(-kappa * r) * (kappa * r + 1) - 1.) / 4.0 / M_PI / r / d;
    };

    auto K1z = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        if (r < vacant_radius) return 0.;
        else return (a.z - b.z) * (dE / dI * exp(-kappa * r) * (kappa * r + 1) - 1.) / 4.0 / M_PI / r / d;
    };

    auto K2 = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);
        if (r < vacant_radius) {
            if (kappa == 0.) return 0.;
            return (exp(-kappa * vacant_radius) - 1.0 + kappa * vacant_radius) / kappa / 2.0 / M_PI / vacant_radius /
                   vacant_radius;
        } else {
            return (1.0 - exp(-kappa * r)) / 4.0 / M_PI / r;
        }
    };

    auto K3xx = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        scalar_t df1 = a.x - b.x;
        scalar_t df2 = a.x - b.x;
        scalar_t delta = 1.0;


        scalar_t t = kappa * r;

        if (r < vacant_radius) {
            return 0.;
        } else {
            return (1 - exp(-t) * (t + 1)) * delta / d / r +
                   (exp(-t) * ((t + 3) * t + 3) - 3.0) * df1 * df2 / d / d / r;
        }
    };

    auto K3yy = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        scalar_t df1 = a.y - b.y;
        scalar_t df2 = a.y - b.y;
        scalar_t delta = 1.0;


        scalar_t t = kappa * r;

        if (r < vacant_radius) {
            return 0.;
        } else {
            return (1 - exp(-t) * (t + 1)) * delta / d / r +
                   (exp(-t) * ((t + 3) * t + 3) - 3.0) * df1 * df2 / d / d / r;
        }
    };

    auto K3zz = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        scalar_t df1 = a.z - b.z;
        scalar_t df2 = a.z - b.z;
        scalar_t delta = 1.0;


        scalar_t t = kappa * r;

        if (r < vacant_radius) {
            return 0.;
        } else {
            return (1 - exp(-t) * (t + 1)) * delta / d / r +
                   (exp(-t) * ((t + 3) * t + 3) - 3.0) * df1 * df2 / d / d / r;
        }
    };

    auto K3xy = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        scalar_t df1 = a.x - b.x;
        scalar_t df2 = a.y - b.y;
        scalar_t delta = 1.0;


        scalar_t t = kappa * r;

        if (r < vacant_radius) {
            return 0.;
        } else {
            return (exp(-t) * ((t + 3) * t + 3) - 3.0) * df1 * df2 / d / d / r;
        }
    };

    auto K3yz = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        scalar_t df1 = a.y - b.y;
        scalar_t df2 = a.z - b.z;
        scalar_t delta = 1.0;


        scalar_t t = kappa * r;

        if (r < vacant_radius) {
            return 0.;
        } else {
            return (exp(-t) * ((t + 3) * t + 3) - 3.0) * df1 * df2 / d / d / r;
        }
    };

    auto K3zx = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        scalar_t df1 = a.z - b.z;
        scalar_t df2 = a.x - b.x;
        scalar_t delta = 1.0;


        scalar_t t = kappa * r;

        if (r < vacant_radius) {
            return 0.;
        } else {
            return (exp(-t) * ((t + 3) * t + 3) - 3.0) * df1 * df2 / d / d / r;
        }
    };


    auto K3yx = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        scalar_t df1 = a.x - b.x;
        scalar_t df2 = a.y - b.y;
        scalar_t delta = 1.0;


        scalar_t t = kappa * r;

        if (r < vacant_radius) {
            return 0.;
        } else {
            return (exp(-t) * ((t + 3) * t + 3) - 3.0) * df1 * df2 / d / d / r;
        }
    };

    auto K3zy = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        scalar_t df1 = a.y - b.y;
        scalar_t df2 = a.z - b.z;
        scalar_t delta = 1.0;


        scalar_t t = kappa * r;

        if (r < vacant_radius) {
            return 0.;
        } else {
            return (exp(-t) * ((t + 3) * t + 3) - 3.0) * df1 * df2 / d / d / r;
        }
    };

    auto K3xz = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        scalar_t df1 = a.z - b.z;
        scalar_t df2 = a.x - b.x;
        scalar_t delta = 1.0;


        scalar_t t = kappa * r;

        if (r < vacant_radius) {
            return 0.;
        } else {
            return (exp(-t) * ((t + 3) * t + 3) - 3.0) * df1 * df2 / d / d / r;
        }
    };


    auto K4x = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        if (r < vacant_radius) return 0.;
        else return -(a.x - b.x) * (dI / dE * exp(-kappa * r) * (kappa * r + 1) - 1.) / 4.0 / M_PI / r / d;
    };

    auto K4y = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        if (r < vacant_radius) return 0.;
        else return -(a.y - b.y) * (dI / dE * exp(-kappa * r) * (kappa * r + 1) - 1.) / 4.0 / M_PI / r / d;
    };

    auto K4z = [&](point &a, point &b) {
        scalar_t d = SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z);
        scalar_t r = sqrt(d);

        if (r < vacant_radius) return 0.;
        else return -(a.z - b.z) * (dI / dE * exp(-kappa * r) * (kappa * r + 1) - 1.) / 4.0 / M_PI / r / d;
    };


    auto dmapping = [&](vector<point> &_source, vector<point> &_target, vector<scalar_t> &_weight,
                        vector<scalar_t> &_normalX, vector<scalar_t> &_normalY, vector<scalar_t> &_normalZ,
                        Vector &_phi) {
        index_t _N = (index_t) _source.size();
        assert(_phi.row() == 2 * _N);
        kernel _K1x, _K1y, _K1z, _K2, _K3xx, _K3yy, _K3zz, _K3xy, _K3yz, _K3zx, _K3yx, _K3xz, _K3zy, _K4x, _K4y, _K4z;

        _K1x.eval = K1x;
        _K1y.eval = K1y;
        _K1z.eval = K1z;

        _K2.eval = K2;
        _K3xx.eval = K3xx;
        _K3yy.eval = K3yy;
        _K3zz.eval = K3zz;
        _K3xy.eval = K3xy;
        _K3yz.eval = K3yz;
        _K3zx.eval = K3zx;
        _K3yx.eval = K3yx;
        _K3zy.eval = K3zy;
        _K3xz.eval = K3xz;

        _K4x.eval = K4x;
        _K4y.eval = K4y;
        _K4z.eval = K4z;


        Vector pphi_pn(_N);
        Vector phiX(_N), phiY(_N), phiZ(_N);

        for (auto id = 0; id < _N; ++id) {
            pphi_pn(id) = _phi(id + _N) * _weight[id];
            phiX(id) = _phi(id) * _normalX[id] * _weight[id];
            phiY(id) = _phi(id) * _normalY[id] * _weight[id];
            phiZ(id) = _phi(id) * _normalZ[id] * _weight[id];
        }

        // todo: check all charges.

        _K1x.initialize(np, _source, _target, phiX, _N, _N, maxPoint, maxLevel);
        _K1y.initialize(np, _source, _target, phiY, _N, _N, maxPoint, maxLevel);
        _K1z.initialize(np, _source, _target, phiZ, _N, _N, maxPoint, maxLevel);
        _K2.initialize(np, _source, _target, pphi_pn, _N, _N, maxPoint, maxLevel);
        _K3xx.initialize(np, _source, _target, phiX, _N, _N, maxPoint, maxLevel);
        _K3yy.initialize(np, _source, _target, phiY, _N, _N, maxPoint, maxLevel);
        _K3zz.initialize(np, _source, _target, phiZ, _N, _N, maxPoint, maxLevel);
        _K3xy.initialize(np, _source, _target, phiX, _N, _N, maxPoint, maxLevel);
        _K3yz.initialize(np, _source, _target, phiY, _N, _N, maxPoint, maxLevel);
        _K3zx.initialize(np, _source, _target, phiZ, _N, _N, maxPoint, maxLevel);
        _K3yx.initialize(np, _source, _target, phiY, _N, _N, maxPoint, maxLevel);
        _K3zy.initialize(np, _source, _target, phiZ, _N, _N, maxPoint, maxLevel);
        _K3xz.initialize(np, _source, _target, phiX, _N, _N, maxPoint, maxLevel);
        _K4x.initialize(np, _source, _target, pphi_pn, _N, _N, maxPoint, maxLevel);
        _K4y.initialize(np, _source, _target, pphi_pn, _N, _N, maxPoint, maxLevel);
        _K4z.initialize(np, _source, _target, pphi_pn, _N, _N, maxPoint, maxLevel);

        Vector ret1x, ret1y, ret1z, ret2, ret3xx, ret3yy, ret3zz, ret3xy, ret3yz, ret3zx, ret3yx, ret3zy, ret3xz, ret4x, ret4y, ret4z;

        _K1x.run(ret1x);
        _K1y.run(ret1y);
        _K1z.run(ret1z);
        _K2.run(ret2);
        _K3xx.run(ret3xx);
        _K3yy.run(ret3yy);
        _K3zz.run(ret3zz);
        _K3xy.run(ret3xy);
        _K3yx.run(ret3yx);
        _K3xz.run(ret3xz);
        _K3zx.run(ret3zx);
        _K3yz.run(ret3yz);
        _K3zy.run(ret3zy);
        _K4x.run(ret4x);
        _K4y.run(ret4y);
        _K4z.run(ret4z);

        Vector tmp1(_N), tmp2(_N);

        for (auto id = 0; id < _N; ++id) {
            tmp1(id) = normalX[id] * (ret3xx(id) + ret3yx(id) + ret3zx(id)) +
                       normalY[id] * (ret3yy(id) + ret3xy(id) + ret3zy(id)) +
                       normalZ[id] * (ret3zz(id) + ret3xz(id) + ret3yz(id));

            tmp2(id) = normalX[id] * ret4x(id) + normalY[id] * ret4y(id) + normalZ[id] * ret4z(id);
        }

        Vector output(2 * _N);

        for (auto id = 0; id < _N; ++id) {
            output(id) = 0.5 * (1 + dE / dI) * _phi(id) + ret1x(id) + ret1y(id) + ret1z(id) - ret2(id);
            output(id + _N) = 0.5 * (1 + dI / dE) * _phi(id + _N) + tmp1(id) - tmp2(id);
        }

        return output;

    };


    auto mapping = [&](vector<point> &_source, vector<point> &_target, vector<scalar_t> &_weight,
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
        return dmapping(source, target, weight, normalX, normalY, normalZ, phi);
    };


    Vector start(2 * N); setValue(start, 0.);

    Vector load(2 * N); setValue(load, 0.);
    for (auto id = 0; id < N; ++id) {
        for (auto atom_id = 0; atom_id <  mol.N; ++atom_id) {
            scalar_t d = SQR(source[id].x - mol.centers[atom_id].data[0] / rescale) +
                         SQR(source[id].y - mol.centers[atom_id].data[1] / rescale) +
                         SQR(source[id].z - mol.centers[atom_id].data[2] / rescale);

            scalar_t r = sqrt(d);
            load(id) += mol.charges[atom_id] / dI / 4.0 / M_PI / r;
            load(id + N) -= mol.charges[atom_id] / dI / 4.0 / M_PI / d / r *
                            (normalX[id] * (source[id].x - mol.centers[atom_id].data[0] / rescale) +
                             normalY[id] * (source[id].y - mol.centers[atom_id].data[1] / rescale) +
                             normalZ[id] * (source[id].z - mol.centers[atom_id].data[2] / rescale));
        }
    }


    GMRES(FullMap, start, load, atoi(cfg.options["gmres_restart"].c_str()), atoi(cfg.options["gmres_maxiter"].c_str()), atof(cfg.options["gmres_tol"].c_str()));

    /*
     * solution is stored in start.
     * output to file.
     */
    std::ofstream potentialFile;
    potentialFile.open(cfg.options["potent_file"]);

    for (int id = 0; id < 2 * N; ++id) {
        potentialFile << start(id) << "\n";
    }
    potentialFile.close();

    /*
     * calculate polarization energy
     */
    vector<point> target_centers;

    for (int id = 0; id < mol.N; ++id) {
        target_centers.push_back(
                {
                    mol.centers[id].data[0]/rescale,
                    mol.centers[id].data[1]/rescale,
                    mol.centers[id].data[2]/rescale
                }
        );

    }

    auto polarizeMap = [&](vector<point>& _source, vector<point>& _target, vector<scalar_t>& _weight,
                           vector<scalar_t>& _normalX,  vector<scalar_t>& _normalY, vector<scalar_t>& _normalZ, Vector& _phi) {

        index_t _N = (index_t) _source.size();
        index_t _M = (index_t) _target.size();

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

        G0.initialize(np, _source, _target, pphi_pn, _N, _M, maxPoint, maxLevel);
        pG0x.initialize(np, _source, _target, phiX, _N, _M, maxPoint, maxLevel);
        pG0y.initialize(np, _source, _target, phiY, _N, _M, maxPoint, maxLevel);
        pG0z.initialize(np, _source, _target, phiZ, _N, _M, maxPoint, maxLevel);

        Gk.initialize(np, _source, _target, pphi_pn, _N, _M, maxPoint, maxLevel);
        pGkx.initialize(np, _source, _target, phiX, _N, _M, maxPoint, maxLevel);
        pGky.initialize(np, _source, _target, phiY, _N, _M, maxPoint, maxLevel);
        pGkz.initialize(np, _source, _target, phiZ, _N, _M, maxPoint, maxLevel);

        Vector retG0, retGk, retG0X, retG0Y, retG0Z, retGkX, retGkY, retGkZ;
        G0.run(retG0); Gk.run(retGk);
        pG0x.run(retG0X); pG0y.run(retG0Y);
        pG0z.run(retG0Z); pGkx.run(retGkX);
        pGky.run(retGkY); pGkz.run(retGkZ);

        Vector output(_M); setValue(output, 0.);

        for (auto id = 0; id < _M; ++id) {
            output(id) = dE / dI * (retGkX(id) + retGkY(id) + retGkZ(id)) - retG0X(id) - retG0Y(id)
                                    - retG0Z(id) + retG0(id) - retGk(id);
        }
        return output;
    };


    Vector energy(mol.N); setValue(energy, 0.);
    energy = polarizeMap(source, target_centers, weight, normalX, normalY, normalZ, start);

    std::ofstream energyFile;
    energyFile.open(cfg.options["energy_file"]);

    for (int id = 0; id < energy.row(); ++id) {
        energyFile << energy(id) << "\n";
    }
    energyFile.close();

    scalar_t polarizedEnergy = 0.;
    for (auto id = 0; id < energy.row(); ++id) {
        polarizedEnergy += energy(id) * mol.charges[id];
    }
    polarizedEnergy *= 0.5;

    std::cout << "polarized energy: " << std::setw(20) << std::scientific <<polarizedEnergy <<std::fixed << std::endl;

}


