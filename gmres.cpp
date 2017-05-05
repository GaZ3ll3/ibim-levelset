//
// Created by lurker on 3/23/17.
//

#include "gmres.h"

using namespace bbfmm;

inline void
Update(Vector &x, int k, Matrix &h, Vector &s, Vector v[]) {
    Vector y(s);

    for (int i = k; i >= 0; i--) {
        y(i) /= h(i, i);
        for (int j = i - 1; j >= 0; j--)
            y(j) -= h(j, i) * y(i);
    }

    for (int j = 0; j <= k; j++)
        daxpy(y(j), v[j], x);
}

void GeneratePlaneRotation(scalar_t &dx, scalar_t &dy, scalar_t &cs, scalar_t &sn) {
    if (dy == 0.0) {
        cs = 1.0;
        sn = 0.0;
    } else if (abs(dy) > abs(dx)) {
        scalar_t temp = dx / dy;
        sn = 1.0 / sqrt(1.0 + temp * temp);
        cs = temp * sn;
    } else {
        scalar_t temp = dy / dx;
        cs = 1.0 / sqrt(1.0 + temp * temp);
        sn = temp * cs;
    }
}


void ApplyPlaneRotation(scalar_t &dx, scalar_t &dy, scalar_t &cs, scalar_t &sn) {
    scalar_t temp = cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
}


inline scalar_t abs(scalar_t x) {
    return (x > 0 ? x : -x);
}

int GMRES(const std::function<Vector(Vector &)> A, Vector &x, Vector &b, int m, int max_iter,
          scalar_t tol) {


    /*
     * for performance use
     */
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    scalar_t resid;
    scalar_t _tol = tol;
    int _max_iter = max_iter;

    int i, j = 1, k;
    Vector s(m + 1), cs(m + 1), sn(m + 1), w;

    scalar_t normb = nrm2(b);
    Vector r = b;
    Vector p = A(x);

    daxpy(-1.0, p, r);

    scalar_t beta = nrm2(r);

    /*
     * error is |Ax - b|/|b|
     */

    std::cout << "=============== GMRES =================" << std::endl;
    std::cout << "    iter    |  rel error   |   time   " << std::endl;

    begin = std::chrono::steady_clock::now();

    if (normb == 0.0)
        normb = 1;

    if ((resid = nrm2(r) / normb) <= _tol) {
        _tol = resid;
        _max_iter = 0;
        return 0;
    }

    Vector *v = new Vector[m + 1];
    Matrix H(x.row(), m + 1);

    while (j <= _max_iter) {

        v[0] = r;    // ??? r / beta
        dscal((1.0 / beta), v[0]);
        setValue(s, 0.);
        s(0) = beta;

        for (i = 0; i < m && j <= _max_iter; i++, j++) {

            end = std::chrono::steady_clock::now();
            std::cout << std::setw(6) << j << std::setw(20) << std::scientific << resid
                      << std::setw(12) << std::fixed
                      << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0 << std::fixed
                      << std::endl;
            begin = std::chrono::steady_clock::now();
            w = A(v[i]);
            for (k = 0; k <= i; k++) {
                H(k, i) = ddot(w, v[k]);
                daxpy(-H(k, i), v[k], w);
            }
            H(i + 1, i) = nrm2(w);
            v[i + 1] = w;
            dscal((1.0 / H(i + 1, i)), v[i + 1]);
            // ??? w / H(i+1, i)

            for (k = 0; k < i; k++)
                ApplyPlaneRotation(H(k, i), H(k + 1, i), cs(k), sn(k));

            GeneratePlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
            ApplyPlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
            ApplyPlaneRotation(s(i), s(i + 1), cs(i), sn(i));

            if ((resid = abs(s(i + 1)) / normb) < _tol) {
                Update(x, i, H, s, v);
                _tol = resid;
                _max_iter = j;

                std::cout << "=============== CONVERGED =============" << std::endl;

                delete[] v;
                return 0;
            }
        }
        Update(x, i - 1, H, s, v);

        r = b;
        p = A(x);

        daxpy(-1.0, p, r);

        beta = nrm2(r);
        if ((resid = beta / normb) < _tol) {
            _tol = resid;
            _max_iter = j;

            std::cout << "=============== CONVERGED =============" << std::endl;

            delete[] v;
            return 0;
        }
    }

    _tol = resid;
    delete[] v;
    std::cout << "============= NOT CONVERGED ===========" << std::endl;
    return 1;
}
