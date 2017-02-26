//
// Created by lurker on 2/23/17.
//

#include "levelset.h"

#define CHUNK 80

Grid::Grid(index_t _Nx, index_t _Ny, index_t _Nz) {
    Nx = _Nx; Ny = _Ny; Nz = _Nz;
    data = new scalar_t[_Nx * _Ny * _Nz];
}

Grid::Grid(scalar_t val, index_t _Nx, index_t _Ny, index_t _Nz) {
    Nx = _Nx; Ny = _Ny; Nz = _Nz;
    data = new scalar_t[_Nx * _Ny * _Nz];
    std::fill_n(data, Nx * Ny * Nz, val);
}

Grid::~Grid() {
#ifdef DEBUG
    printf("deconstructing grid...\n");
#endif
    delete[] data;
    data = nullptr;
}

scalar_t Grid::get(index_t i, index_t j, index_t k) {
    i = min(max(i, 0), Nx - 1);
    j = min(max(j, 0), Ny - 1);
    k = min(max(k, 0), Nz - 1);
    return data[i * Ny * Nz + j * Nz + k];
}

void Grid::set(scalar_t val, index_t i, index_t j, index_t k) {
    data[i * Ny * Nz + j * Nz + k] = val;
}



levelset::levelset(index_t _Nx, index_t _Ny, index_t _Nz, index_t _band, scalar_t _sx, scalar_t _sy, scalar_t _sz, scalar_t _dx) {
    dx = _dx; Nx = _Nx; Ny = _Ny; Nz = _Nz;
    sx = _sx; sy = _sy; sz = _sz; bandwidth = _band;
    shift = 7; // for WENO5.
}


levelset::~levelset() {}

void levelset::expand(Molecule &mol, Grid &g, scalar_t probe) {
    for (index_t aId = 0; aId < (index_t)mol.centers.size();++aId) {
        int r = (int)((mol.radii[aId] + probe)/dx) + bandwidth;
        int icx = (int)((mol.centers[aId].data[0] - sx) / dx);
        int icy = (int)((mol.centers[aId].data[1] - sy) / dx);
        int icz = (int)((mol.centers[aId].data[2] - sz) / dx);
#pragma omp parallel for schedule(static) collapse(3) num_threads(4)
        for (int a = icx - r; a <= icx + r; ++a) {
            for (int b = icy - r; b <= icy + r; ++b) {
                for (int c = icz - r; c<= icz + r; ++c) {
                    point grid_p = {sx + a * dx, sy + b * dx, sz + c * dx};
                    scalar_t dist = (mol.radii[aId] + probe) - norm(grid_p - mol.centers[aId]) ;
                    dist = max(dist, g.get(a, b,c));
                    g.set(dist, a, b, c);
                }
            }
        }
    }
}

void levelset::evolve(Grid &g, scalar_t final_t, scalar_t vel, scalar_t cfl_thres) {
    scalar_t dt = cfl_thres * dx / vel;
    index_t num_steps = (index_t)(final_t / dt) + 1;
    dt = final_t / (scalar_t)num_steps;

    std::cout << std::setw(15) << "INWARD CONFIG" << std::endl;

    std::cout << std::setw(15)<< "dt" << " " << std::setw(8) << dt << " seconds" <<std::endl;
    std::cout << std::setw(15)<< "STEPS" << " " << std::setw(8) << num_steps << " steps" <<std::endl;


    Grid u1(Nx, Ny, Nz);
    Grid u2(Nx, Ny, Nz);

    index_t step = 0;
    point Dup, Dun;

    auto core = omp_get_max_threads();
    omp_set_num_threads(core);

    double* window = (double*)malloc((core)*DIM * shift * sizeof(double));

    for (step = 0; step < num_steps; ++step) {
#pragma omp parallel for private(Dup, Dun) schedule(static, CHUNK) collapse(3) num_threads(core)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;
                    index_t tid = omp_get_thread_num();

                    setWindow(g, window + tid * DIM * shift, i, j, k);

                    for (index_t dir = 0; dir < DIM; ++dir) {
                        setGradient(dir, window + tid * DIM * shift, Dup, Dun);
                    }

                    u1.data[I] = g.data[I] - dt * vel * getNorm(Dun, Dup);
                }
            }
        }
#pragma omp barrier
#pragma omp parallel for private(Dup, Dun) schedule(static, CHUNK) collapse(3) num_threads(core)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;
                    index_t tid = omp_get_thread_num();

                    setWindow(u1, window+ tid * DIM * shift, i, j, k);

                    for (index_t dir = 0; dir < DIM; ++dir) {
                        setGradient(dir,window+ tid * DIM * shift, Dup, Dun);
                    }

                    u2.data[I] = (3 * g.data[I] + u1.data[I] - dt * vel * getNorm(Dun, Dup)) / 4.0;
                }
            }
        }
#pragma omp barrier
#pragma omp parallel for private(Dup, Dun) schedule(static, CHUNK) collapse(3) num_threads(core)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;
                    index_t tid = omp_get_thread_num();

                    setWindow(u2, window+ tid * DIM * shift, i, j, k);

                    for (index_t dir = 0; dir < DIM; ++dir) {
                        setGradient(dir, window+ tid * DIM * shift, Dup, Dun);
                    }

                    g.data[I] =( g.data[I] + 2 * (u2.data[I] - dt * vel * getNorm(Dun, Dup)) ) / 3.0;
                }
            }
        }
#pragma omp barrier
    }
    free(window);
}



void levelset::reinitialize(Grid &g, Grid &phi0, scalar_t final_t, scalar_t vel, scalar_t cfl_thres) {
    scalar_t dt = cfl_thres * dx / vel;
    scalar_t eps = dx * dx;
    index_t num_steps = (index_t)(final_t / dt) + 1;
    dt = final_t / (scalar_t)num_steps;

    std::cout << std::setw(15) << "REINIT CONFIG" << std::endl;

    std::cout << std::setw(15)<< "dt" << " " << std::setw(8) << dt << " seconds" <<std::endl;
    std::cout << std::setw(15)<< "STEPS" << " " << std::setw(8) << num_steps << " steps" <<std::endl;


    Grid u1(Nx, Ny, Nz);
    Grid u2(Nx, Ny, Nz);

    index_t step = 0;
    point Dup, Dun;

    auto core = omp_get_max_threads();
    omp_set_num_threads(core);

    double* window = (double*)malloc((core)*DIM * shift * sizeof(double));

    int indices = 0;
    int total = 0;
    int thickness = 3;
    scalar_t thres = 1e-2;

    scalar_t *_window = (scalar_t *)malloc(DIM * shift * sizeof(scalar_t));
    for (index_t i = 0; i < Nx; ++i) {
        for (index_t j = 0; j < Ny; ++j) {
            for (index_t k= 0; k < Nz; ++k) {

                if ( fabs(g.get(i, j, k)) < thickness * dx ) {
                    total++;
                    point _Dun, _Dup;
                    setWindow(g, _window, i, j, k);
                    for (index_t dir = 0; dir < DIM; ++dir) {
                        setGradient(dir, _window, _Dup, _Dun);
                    }

                    scalar_t nr = getNorm(_Dun, _Dup);

                    if (fabs(nr - 1.0) > thres) {
                        indices++;
                    }
                }
            }
        }
    }
    std::cout << "incorrect gradients :" << indices << " , " << total << std::endl;

    while (indices !=0 && step < num_steps) {
        step++;
#pragma omp parallel for private(Dup, Dun) schedule(static, CHUNK) collapse(3) num_threads(core)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;
                    index_t tid = omp_get_thread_num();

                    setWindow(g, window + tid * DIM * shift, i, j, k);

                    for (index_t dir = 0; dir < DIM; ++dir) {
                        setGradient(dir, window + tid * DIM * shift, Dup, Dun);
                    }

                    double sign = phi0.data[I] / (sqrt(SQR(phi0.data[I]) + eps));
                    double normDu = (sign > 0. ? getNorm(Dun, Dup) : getNorm(Dup, Dun));

                    u1.data[I] = g.data[I] - dt * vel * sign * (normDu - 1.0);
                }
            }
        }
#pragma omp barrier
#pragma omp parallel for private(Dup, Dun) schedule(static, CHUNK) collapse(3) num_threads(core)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;
                    index_t tid = omp_get_thread_num();

                    setWindow(u1, window+ tid * DIM * shift, i, j, k);

                    for (index_t dir = 0; dir < DIM; ++dir) {
                        setGradient(dir,window+ tid * DIM * shift, Dup, Dun);
                    }

                    double sign = phi0.data[I] / (sqrt(SQR(phi0.data[I]) + eps));
                    double normDu = (sign > 0. ? getNorm(Dun, Dup) : getNorm(Dup, Dun));

                    u2.data[I] = (3 * g.data[I] + u1.data[I] - dt * vel * sign *  (normDu - 1.0)) / 4.0;
                }
            }
        }
#pragma omp barrier
#pragma omp parallel for private(Dup, Dun) schedule(static, CHUNK) collapse(3) num_threads(core)
        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k = 0; k < Nz; ++k) {
                    index_t I = i * Ny * Nz + j * Nz + k;
                    index_t tid = omp_get_thread_num();

                    setWindow(u2, window+ tid * DIM * shift, i, j, k);

                    for (index_t dir = 0; dir < DIM; ++dir) {
                        setGradient(dir, window+ tid * DIM * shift, Dup, Dun);
                    }

                    double sign = phi0.data[I] / (sqrt(SQR(phi0.data[I]) + eps));
                    double normDu = (sign > 0. ? getNorm(Dun, Dup) : getNorm(Dup, Dun));

                    g.data[I] =( g.data[I] + 2 * (u2.data[I] - dt * vel * sign *  (normDu - 1.0)) ) / 3.0;
                }
            }
        }
#pragma omp barrier

        scalar_t *_window = (scalar_t *)malloc(DIM * shift * sizeof(scalar_t));

        indices = 0;
        total = 0;

        for (index_t i = 0; i < Nx; ++i) {
            for (index_t j = 0; j < Ny; ++j) {
                for (index_t k= 0; k < Nz; ++k) {

                    if ( fabs(g.get(i, j, k)) < thickness * dx) {
                        total++;
                        point _Dun, _Dup;
                        setWindow(g, _window, i, j, k);
                        for (index_t dir = 0; dir < DIM; ++dir) {
                            setGradient(dir, _window, _Dup, _Dun);
                        }

                        scalar_t nr = getNorm(_Dun, _Dup);

                        if (fabs(nr - 1.0) > thres) {
                            indices++;
                        }
                    }
                }
            }
        }
        std::cout << "incorrect gradients :" << indices << " , " << total << std::endl;
        free(_window);


    }
    free(window);

    vector<index_t > IND;

    total = 0;

    for (index_t i = 0; i < Nx; ++i) {
        for (index_t j = 0; j < Ny; ++j) {
            for (index_t k= 0; k < Nz; ++k) {

                if ( fabs(g.get(i, j, k)) < thickness * dx ) {
                    total++;
                    point _Dun, _Dup;
                    setWindow(g, _window, i, j, k);
                    for (index_t dir = 0; dir < DIM; ++dir) {
                        setGradient(dir, _window, _Dup, _Dun);
                    }

                    scalar_t nr = getNorm(_Dun, _Dup);

                    if (fabs(nr - 1.0) > thres) {
                        IND.push_back(i * Ny * Nz + j * Nz + k);
                    }
                }
            }
        }
    }


    for (index_t i = 0; i < Nx; ++i) {
        for (index_t j = 0; j < Ny; ++j) {
            for (index_t k = 0; k < Nz; ++k) {
                g.data[i * Ny * Nz + j* Nz + k] = 10000.0;
            }
        }
    }

    for (auto I : IND) {
        g.data[I] = 0.;
    }

    free(_window);


}



scalar_t levelset::getNorm(point &Dun, point &Dup) {
    scalar_t val = 0.;
    for (int i = 0; i < DIM; ++i) {
        Dun.data[i] = max((Dun.data[i] > 0. ? Dun.data[i] : 0.), (Dup.data[i] < 0. ? -Dup.data[i] : 0.));
        val += SQR(Dun.data[i]);
    }
    return sqrt(val);
}

void levelset::setWindow(Grid &g, scalar_t *window, index_t i, index_t j, index_t k) {
    index_t half = shift / 2;
    index_t _loc = 0;
    index_t xs = i - half, ys = j - half, zs = k - half;
    for (index_t loc = 0; loc < shift; ++loc) {
        window[_loc++] = g.get(xs + loc, j, k);
    }
    for (index_t loc = 0; loc < shift; ++loc) {
        window[_loc++] = g.get(i, ys+loc, k);
    }

    for (index_t loc = 0; loc < shift; ++loc) {
        window[_loc++] = g.get(i, j, zs+loc);
    }
}

void levelset::setGradient(index_t dir, scalar_t *window, point &uxp, point &uxn) {
    scalar_t df, is0, is1, is2, a0, a1, a2, w0, w1, w2;
    scalar_t a,b,c,d,e,f;
    scalar_t EPS_LS = 1E-6;

    index_t start = dir * shift;

    a = *(window + start + 1) - *(window + start);
    b = *(window + start + 2) - *(window + start + 1);
    c = *(window + start + 3) - *(window + start + 2);
    d = *(window + start + 4) - *(window + start + 3);
    e = *(window + start + 5) - *(window + start + 4);
    f = *(window + start + 6) - *(window + start + 5);


    df = (-b + 7.0 * (c + d) - e) / 12.0;

    a = b - a;
    b = c - b;
    c = d - c;
    d = e - d;
    e = f - e;

    is0 = 13.0 * SQR(e - d) + 3.0 * SQR(e- 3.0 * d);
    is1 = 13.0 * SQR(d - c) + 3.0 * SQR(d + c);
    is2 = 13.0 * SQR(c - b) + 3.0 * SQR(3.0 * c- b);

    a0 = 1.0 / SQR(EPS_LS + is0);
    a1 = 6.0 / SQR(EPS_LS + is1);
    a2 = 3.0 / SQR(EPS_LS + is2);

    w1 = 1.0 / (a0 + a1 + a2);
    w0 = a0 * w1;
    w2 = a2 * w1;

    uxp.data[dir] = df +
                    w0 * (e - 2 * d + c) / 3.0 +
                    (w2 - 0.5) * (d - 2.0 * c + b) / 6.0;

    uxp.data[dir] /= dx;

    is0 = 13.0 * SQR(a - b) + 3.0 * SQR(a - 3.0 * b);
    is1 = 13.0 * SQR(b - c) + 3.0 * SQR(b + c);
    is2 = 13.0 * SQR(c - d) + 3.0 * SQR(3.0 * c - d);

    a0 = 1.0 / SQR(EPS_LS + is0);
    a1 = 6.0 / SQR(EPS_LS + is1);
    a2 = 3.0 / SQR(EPS_LS + is2);

    w1 = 1.0 / (a0 + a1 + a2);
    w0 = a0 * w1;
    w2 = a2 * w1;

    uxn.data[dir] = df -
                    w0 * (a - 2 * b + c) / 3.0 -
                    (w2 - 0.5) * (b - 2 * c + d) / 6.0;

    uxn.data[dir] /= dx;
}


