//
// Created by lurker on 2/10/17.
//

#ifndef POINT_H
#define POINT_H

#include "utils.h"

class ls_basePoint {
public:
    scalar_t data[DIM];

    ls_basePoint()  {memset(data, 0, DIM * sizeof(scalar_t));}

    ls_basePoint(scalar_t _x, scalar_t _y, scalar_t _z) {
        data[0] = _x;
        data[1] = _y;
        data[2] = _z;
    }

    virtual ~ls_basePoint() {}

    bool operator>=(const ls_basePoint &a) {
        return (data[0] >= a.data[0] - EPS) && (data[1] >= a.data[1] - EPS) && (data[2] >= a.data[2] - EPS);
    }

    bool operator<=(const ls_basePoint &a) {
        return (data[0] <= a.data[0] + EPS) && (data[1] <= a.data[1] + EPS) && (data[2] <= a.data[2] + EPS);
    }

    bool operator==(const ls_basePoint &a) {
        return fabs(data[0] - a.data[0]) < EPS && fabs(data[1] - a.data[1]) < EPS && fabs(data[2] - a.data[2]) < EPS;
    }
};

class ls_point : public ls_basePoint {
public:
    ls_point() : ls_basePoint() {
    }

    ls_point(scalar_t _x, scalar_t _y, scalar_t _z) : ls_basePoint(_x, _y, _z) {
    }

    ls_point(const ls_point& p) {
        cblas_dcopy(DIM, p.data, 1, this->data, 1);
    }

    ls_point(ls_point&& p) {
        cblas_dcopy(DIM, p.data, 1, this->data, 1);
    }

    ~ls_point() {}

    /*
     * linear algebra for vec3.
     */
    ls_point operator+ (const ls_point& p) {
        ls_point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, 1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    ls_point operator+(const scalar_t val) {
        ls_point p; std::fill_n(p.data, DIM, val);
        ls_point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, 1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    ls_point operator-(const ls_point& p) {
        ls_point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, -1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    ls_point operator+ (ls_point&& p) {
        ls_point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, 1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    ls_point operator-(const scalar_t val) {
        ls_point p; std::fill_n(p.data, DIM, val);
        ls_point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, -1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    ls_point operator-(ls_point&& p) {
        ls_point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, -1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    ls_point operator*(scalar_t val) {
        ls_point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_dscal(DIM, val, ret.data, 1);
        return ret;
    }

    void operator=(const ls_point& p) {
        cblas_dcopy(DIM, p.data, 1, this->data, 1);
    }

    void operator=(ls_point&& p) {
        cblas_dcopy(DIM, p.data, 1, this->data, 1);
    }
};

inline std::ostream &operator<<(std::ostream &os, ls_point &p) {
    for (int i = 0; i < DIM - 1; ++i) {
        os << p.data[i] << ", ";
    }
    os << p.data[DIM - 1];
    return os;
}

inline scalar_t norm(const ls_point& a) {
    return sqrt(cblas_ddot(DIM, a.data, 1, a.data, 1));
}

inline scalar_t nrm2(const ls_point& a) {
    return cblas_ddot(DIM, a.data, 1, a.data, 1);
}

inline scalar_t norm(ls_point&& a) {
    return sqrt(cblas_ddot(DIM, a.data, 1, a.data, 1));
}

inline scalar_t nrm2(ls_point&& a) {
    return cblas_ddot(DIM, a.data, 1, a.data, 1);
}

#endif //POINT_H
