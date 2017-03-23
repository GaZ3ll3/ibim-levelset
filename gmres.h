//
// Created by lurker on 3/22/17.
//

#ifndef LEVELSET_GMRES_H
#define LEVELSET_GMRES_H

#include "blas_wraper.h"

inline void
Update(Vector &x, int k, Matrix &h, Vector &s, Vector v[]) ;

void GeneratePlaneRotation(scalar_t &dx, scalar_t &dy, scalar_t &cs, scalar_t &sn) ;

void ApplyPlaneRotation(scalar_t &dx, scalar_t &dy, scalar_t &cs, scalar_t &sn) ;


inline scalar_t abs(scalar_t x);

int GMRES(const std::function<Vector(Vector &)> A, Vector &x, Vector &b, int m, int max_iter,
          scalar_t tol);

#endif //LEVELSET_GMRES_H
