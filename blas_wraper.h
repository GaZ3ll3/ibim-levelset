//
// Created by lurker on 3/22/17.
//

#ifndef LEVELSET_LINALG_AUX_H
#define LEVELSET_LINALG_AUX_H


#include "linalg.h"

using namespace bbfmm;
/*
 * wrapper of blas.
*/
/*
 * x = ax
 */
void dscal(scalar_t alpha, Vector &v) ;

void dscal(scalar_t alpha, Matrix &v) ;

/*
 *  y = ax + y
 */
void daxpy(scalar_t a, Vector &x, Vector &y);

/*
 * unsafe case
 */
void daxpy(scalar_t a, Vector &x, scalar_t *y) ;

/*
 *  y = ax + y
 */
void daxpy(scalar_t a, Matrix &x, Matrix &y) ;

/*
 * C = aAB + bC
 */
void dgemm(scalar_t alpha, Matrix &A, Matrix &B, scalar_t beta, Matrix &C) ;
/*
 * C = aABt + bC
 */
void dgemm_t(scalar_t alpha, Matrix &A, Matrix &B, scalar_t beta, Matrix &C) ;

/*
 * C = aAtB + bC
 */
void t_dgemm(scalar_t alpha, Matrix &A, Matrix &B, scalar_t beta, Matrix &C) ;


/*
 * rank-1 update
 * A = a x * y' + A
 */
void dger(scalar_t alpha, Vector &X, Vector &Y, Matrix &A) ;

/*
 * y = a Ax + b y
 */
void dgemv(scalar_t alpha, Matrix &A, Vector &x, scalar_t beta, Vector &y) ;

/*
 * y = a Atx + b y
 */
void dgemv_t(scalar_t alpha, Matrix &A, Vector &x, scalar_t beta, Vector &y) ;

/*
 * hadamard product
 */
void dsbmv(scalar_t alpha, Vector &A, Vector &x, scalar_t beta, Vector &y) ;

/*
 * unsafe case
 */
void dsbmv(scalar_t alpha, Vector &A, Vector &x, scalar_t beta, scalar_t *y) ;

scalar_t nrm2(Vector &x) ;

scalar_t ddot(Vector &x, Vector &y);

#endif