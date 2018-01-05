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

void Grid::output(std::string filename) {
    std::ofstream gridFile;
    gridFile.open(filename);

    for (index_t i = 0; i < Nx; ++i) {
        for (index_t j = 0; j < Ny; ++j) {
            for (index_t k = 0; k < Nz; ++k) {
                index_t I = i * Ny * Nz + j * Nz + k;
                gridFile << data[I] << "\n";
            }
        }
    }
    gridFile.close();
}

Grid &Grid::operator*=(const scalar_t multiplier) {
    for (int i = 0; i < Nx * Ny * Nz; ++i) {
        this->data[i] *= multiplier;
    }
    return *this;
}

Grid &Grid::operator=(const Grid &other) {
    // check size;
    assert(other.Nx == this->Nx);
    assert(other.Ny == this->Ny);
    assert(other.Nz == this->Nz);

    for (int i = 0; i < Nx * Ny * Nz; ++i) {
        this->data[i] = other.data[i];
    }

    return *this;

}


levelset::levelset(index_t _Nx, index_t _Ny, index_t _Nz, index_t _band, scalar_t _sx, scalar_t _sy, scalar_t _sz, scalar_t _dx, Config& cfg) {
    dx = _dx; Nx = _Nx; Ny = _Ny; Nz = _Nz;
    sx = _sx; sy = _sy; sz = _sz; bandwidth = _band;
    shift = 7; // for WENO5.
    thickness = atoi(cfg.options["lvlset_thickness"].c_str());
    thres = atof(cfg.options["lvlset_thres"].c_str());
}


levelset::~levelset() {
#ifdef DEBUG
    std::cout << "destructing levelset." << std::endl;
#endif
}

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
                    ls_point grid_p = {sx + a * dx, sy + b * dx, sz + c * dx};
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

//    std::cout << std::setw(15)<< "INWARD" << " " << std::setw(8) << num_steps << " steps" <<std::endl;


    Grid u1(Nx, Ny, Nz);
    Grid u2(Nx, Ny, Nz);

    index_t step = 0;
    ls_point Dup, Dun;

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
    scalar_t eps = dx * dx; //todo:
    index_t num_steps = (index_t)(final_t / dt) + 1;
    dt = final_t / (scalar_t)num_steps;

//    std::cout << std::setw(15)<< "REINIT" << " " << std::setw(8) << num_steps << " steps" <<std::endl;

    Grid u1(Nx, Ny, Nz);
    Grid u2(Nx, Ny, Nz);

    index_t step = 0;
    ls_point Dup, Dun;

    auto core = omp_get_max_threads();
    omp_set_num_threads(core);

    double* window = (double*)malloc((core)*DIM * shift * sizeof(double));

    int indices = 0;

    indices = countGradient(g, thickness, thres, window);

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
        indices = countGradient(g, thickness, thres, window, false);

    }

    countGradient(g, thickness, thres, window, true);
    free(window);


}


/// getNorm
/// \param Dun : input backward gradient vector.
/// \param Dup : input forward gradient vector.
/// \return : norm
scalar_t levelset::getNorm(ls_point &Dun, ls_point &Dup) {
    scalar_t val = 0.;
    for (int i = 0; i < DIM; ++i) {
        Dun.data[i] = max((Dun.data[i] > 0. ? Dun.data[i] : 0.), (Dup.data[i] < 0. ? -Dup.data[i] : 0.));
        val += SQR(Dun.data[i]);
    }
    return sqrt(val);
}

/// setWindow
/// \param g : input grid
/// \param window : window pointer
/// \param i : input x axis
/// \param j : input y axis
/// \param k : input z axis
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

/// setGradient
/// \param dir : input axis, 0, 1, 2 for x, y, z.
/// \param window : window pointer.
/// \param uxp : output forward gradient
/// \param uxn : output backward gradient
void levelset::setGradient(index_t dir, scalar_t *window, ls_point &uxp, ls_point &uxn) {
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

/// countGradient
/// \param g grid
/// \param thickness tube layer number
/// \param thres gradient pass threshold
/// \param _window shifting window for WENO
/// \param display for displaying information.
/// \return
// todo: remove redundant input.
index_t levelset::countGradient(Grid &g, scalar_t thickness, scalar_t thres, scalar_t* _window, bool display) {
    ls_point _Dun, _Dup;

    index_t total = 0;
    index_t indices = 0;
    scalar_t accum_error = 0.;

    for (index_t i = 0; i < Nx; ++i) {
        for (index_t j = 0; j < Ny; ++j) {
            for (index_t k= 0; k < Nz; ++k) {
                if ( fabs(g.get(i, j, k)) < thickness * dx ) {
                    total++;
                    setWindow(g, _window, i, j, k);
                    for (index_t dir = 0; dir < DIM; ++dir) {
                        setGradient(dir, _window, _Dup, _Dun);
                    }

                    scalar_t nr = getNorm(_Dun, _Dup);

                    if (fabs(nr - 1.0) > thres) {
                            indices++;
                    }
                    accum_error += fabs(nr - 1.0);
                }
            }
        }
    }
    if (display) {
        std::cout << std::setw(15) << "QUALITY" << " " << std::setprecision(2) << std::setw(8)
                  << scalar_t(accum_error) / scalar_t(total) << " & " << scalar_t(indices) / scalar_t(total)
                  << std::endl;
    }

    return indices;
}

/// build surface, remove inclusions.
/// \param g
/// \param ls
Surface::Surface(Grid &g, levelset &ls, scalar_t s) {

    ls_point _Dun, _Dup;

    scalar_t* _window =  (double*)malloc(DIM * ls.shift * sizeof(double));
    scalar_t tube_width = ls.thickness * ls.dx;

    vector<short> visited((unsigned long) (ls.Nx * ls.Ny * ls.Nz), 0);

    queue<int> Queue;
    Queue.push(0);
    visited[0] = 1;

    while (!Queue.empty()) {
        int cur_index = Queue.front();

        Queue.pop();

        int I = cur_index / (ls.Ny * ls.Nz);
        int J = cur_index / ls.Nz - I * ls.Ny;
        int K = cur_index % ls.Nz;
        /*
         * 6 directions to go
         */
        int cI, cJ, cK, next_index;

        cI = I + 1;
        cJ = J;
        cK = K;
        if (cI <= ls.Nx - 1) {
            next_index = cI * ls.Ny * ls.Nz + cJ * ls.Nz + cK;
            if (!visited[next_index] && g.data[cur_index] > -tube_width) {
                visited[next_index] = 1;
                Queue.push(next_index);
            }
        }

        cI = I - 1;
        cJ = J;
        cK = K;
        if (cI >= 0) {
            next_index = cI * ls.Ny * ls.Nz + cJ * ls.Nz + cK;
            if (!visited[next_index] && g.data[cur_index] > -tube_width) {
                visited[next_index] = 1;
                Queue.push(next_index);
            }
        }
        cI = I;
        cJ = J + 1;
        cK = K;
        if (cJ <= ls.Ny - 1) {
            next_index = cI * ls.Ny * ls.Nz + cJ * ls.Nz + cK;
            if (!visited[next_index] && g.data[cur_index] > -tube_width) {
                visited[next_index] = 1;
                Queue.push(next_index);
            }
        }
        cI = I;
        cJ = J - 1;
        cK = K;
        if (cJ >= 0) {
            next_index = cI * ls.Ny * ls.Nz + cJ * ls.Nz + cK;
            if (!visited[next_index] && g.data[cur_index] > -tube_width) {
                visited[next_index] = 1;
                Queue.push(next_index);
            }
        }
        cI = I;
        cJ = J;
        cK = K + 1;
        if (cK <= ls.Nz - 1) {
            next_index = cI * ls.Ny * ls.Nz + cJ * ls.Nz + cK;
            if (!visited[next_index] && g.data[cur_index] > -tube_width) {
                visited[next_index] = 1;
                Queue.push(next_index);
            }
        }
        cI = I;
        cJ = J;
        cK = K - 1;
        if (cK >= 0) {
            next_index = cI * ls.Ny * ls.Nz + cJ * ls.Nz + cK;
            if (!visited[next_index] && g.data[cur_index] > -tube_width) {
                visited[next_index] = 1;
                Queue.push(next_index);
            }
        }
    }

    /*
     * set a value lying outside tube, epsilon
     */
    auto inclusion_value = -(tube_width + 1e-15);
    int cnt = 0;
    for (index_t id = 0; id < ls.Nx * ls.Ny * ls.Nz; ++id) {
        if (!visited[id] && fabs(g.data[id]) < tube_width) {
            g.data[id] = inclusion_value;
            cnt++;
        }
    }
    std::cout << std::setw(15) << "INCLUSION" << " " << std::setw(8) << cnt << std::endl;


    for (index_t i = 0; i < ls.Nx; ++i) {
        for (index_t j = 0; j < ls.Ny; ++j) {
            for (index_t k= 0; k < ls.Nz; ++k) {
                if ( fabs(g.get(i, j, k)) <  tube_width) {
                    ls.setWindow(g, _window, i, j, k);
                    for (index_t dir = 0; dir < DIM; ++dir) {
                        ls.setGradient(dir, _window, _Dup, _Dun);
                    }

                    /*
                     * current gradient is qualified.
                     */
                    _Dun = (_Dun + _Dup) * 0.5;

                    ls_point P = {
                            ls.sx + i * ls.dx,
                            ls.sy + j * ls.dx,
                            ls.sz + k * ls.dx
                    };

                    scalar_t dist = g.get(i, j, k);
                    /*
                     * projection
                     */
                    P = P - _Dun * dist;

                    nodes.push_back(P);
                    normals.push_back(_Dun);
                    /*
                     * calculates the weight according to the distance.
                     *
                     * 1. cos weight.
                     * 2. triangle weight
                     * 3. infinitely smooth weight
                     *
                     */
                    weight.push_back(0.5 * (1.0 + cos(M_PI * dist/tube_width)) / tube_width);
//                    weight.push_back((1.0 - fabs(dist)/tube_width)/tube_width);


                }
            }
        }
    }

    free(_window);
}

Surface::~Surface() {

}


void Surface::output(std::string filename) {
    std::ofstream nodeFile;
    nodeFile.open(filename);

    for (index_t i = 0; i < nodes.size(); ++i) {
        nodeFile << nodes[i].data[0] << " " << nodes[i].data[1] << " " << nodes[i].data[2] << "\n";
    }
    nodeFile.close();
}
