//
// Created by lurker on 2/10/17.
//

#ifndef POINT_H
#define POINT_H

#include "utils.h"

class basePoint {
public:
    scalar_t data[DIM];

    basePoint()  {memset(data, 0, DIM * sizeof(scalar_t));}

    basePoint(scalar_t _x, scalar_t _y, scalar_t _z) {
        data[0] = _x;
        data[1] = _y;
        data[2] = _z;
    }

    virtual ~basePoint() {}

    bool operator>=(const basePoint &a) {
        return (data[0] >= a.data[0] - EPS) && (data[1] >= a.data[1] - EPS) && (data[2] >= a.data[2] - EPS);
    }

    bool operator<=(const basePoint &a) {
        return (data[0] <= a.data[0] + EPS) && (data[1] <= a.data[1] + EPS) && (data[2] <= a.data[2] + EPS);
    }

    bool operator==(const basePoint &a) {
        return fabs(data[0] - a.data[0]) < EPS && fabs(data[1] - a.data[1]) < EPS && fabs(data[2] - a.data[2]) < EPS;
    }
};

class point : public basePoint {
public:
    point() : basePoint() {
    }

    point(scalar_t _x, scalar_t _y, scalar_t _z) : basePoint(_x, _y, _z) {
    }

    point(const point& p) {
        cblas_dcopy(DIM, p.data, 1, this->data, 1);
    }

    point(point&& p) {
        cblas_dcopy(DIM, p.data, 1, this->data, 1);
    }

    ~point() {}

    /*
     * linear algebra for vec3.
     */
    point operator+ (const point& p) {
        point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, 1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    point operator+(const scalar_t val) {
        point p; std::fill_n(p.data, DIM, val);
        point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, 1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    point operator-(const point& p) {
        point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, -1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    point operator+ (point&& p) {
        point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, 1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    point operator-(const scalar_t val) {
        point p; std::fill_n(p.data, DIM, val);
        point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, -1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    point operator-(point&& p) {
        point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_daxpy(DIM, -1.0, p.data, 1, ret.data, 1);
        return ret;
    }

    point operator*(scalar_t val) {
        point ret;
        cblas_dcopy(DIM, this->data, 1, ret.data, 1);
        cblas_dscal(DIM, val, ret.data, 1);
        return ret;
    }

    void operator=(const point& p) {
        cblas_dcopy(DIM, p.data, 1, this->data, 1);
    }

    void operator=(point&& p) {
        cblas_dcopy(DIM, p.data, 1, this->data, 1);
    }
};

inline std::ostream &operator<<(std::ostream &os, point &p) {
    for (int i = 0; i < DIM - 1; ++i) {
        os << p.data[i] << ", ";
    }
    os << p.data[DIM - 1];
    return os;
}

inline scalar_t norm(const point& a) {
    return sqrt(cblas_ddot(DIM, a.data, 1, a.data, 1));
}

inline scalar_t nrm2(const point& a) {
    return cblas_ddot(DIM, a.data, 1, a.data, 1);
}

inline scalar_t norm(point&& a) {
    return sqrt(cblas_ddot(DIM, a.data, 1, a.data, 1));
}

inline scalar_t nrm2(point&& a) {
    return cblas_ddot(DIM, a.data, 1, a.data, 1);
}

#endif //POINT_H
