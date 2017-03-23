//
// Created by lurker on 3/23/17.
//

#include "blas_wraper.h"

using namespace bbfmm;
/*
 * wrapper of blas.
*/
/*
 * x = ax
 */
void dscal(scalar_t alpha, Vector &v) {
    cblas_dscal(v.row(), alpha, v.data(), 1);
}

void dscal(scalar_t alpha, Matrix &v) {
    cblas_dscal(v.row() * v.col(), alpha, v.data(), 1);
}

/*
 *  y = ax + y
 */
void daxpy(scalar_t a, Vector &x, Vector &y) {
    assert(x.row() == y.row());
    cblas_daxpy(x.row(), a, x.data(), 1, y.data(), 1);
}

/*
 * unsafe case
 */
void daxpy(scalar_t a, Vector &x, scalar_t *y) {
    cblas_daxpy(x.row(), a, x.data(), 1, y, 1);
}

/*
 *  y = ax + y
 */
void daxpy(scalar_t a, Matrix &x, Matrix &y) {
    assert(x.row() == y.row());
    assert(x.col() == y.col());
    cblas_daxpy(x.row() * x.col(), a, x.data(), 1, y.data(), 1);
}

/*
 * C = aAB + bC
 */
void dgemm(scalar_t alpha, Matrix &A, Matrix &B, scalar_t beta, Matrix &C) {
    assert(A.row() == C.row());
    assert(A.col() == B.row());
    assert(B.col() == C.col());
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                C.row(), C.col(), A.col(), alpha, A.data(), A.row(), B.data(), B.row(), beta, C.data(), C.row());
}

/*
 * C = aABt + bC
 */
void dgemm_t(scalar_t alpha, Matrix &A, Matrix &B, scalar_t beta, Matrix &C) {
    assert(A.row() == C.row());
    assert(A.col() == B.col());
    assert(B.row() == C.col());
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                C.row(), C.col(), A.col(), alpha, A.data(), A.row(), B.data(), B.row(), beta, C.data(), C.row());
}

/*
 * C = aAtB + bC
 */
void t_dgemm(scalar_t alpha, Matrix &A, Matrix &B, scalar_t beta, Matrix &C) {
    assert(A.col() == C.row());
    assert(A.row() == B.row());
    assert(B.col() == C.col());
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                C.row(), C.col(), A.col(), alpha, A.data(), A.row(), B.data(), B.row(), beta, C.data(), C.row());
}


/*
 * rank-1 update
 * A = a x * y' + A
 */
void dger(scalar_t alpha, Vector &X, Vector &Y, Matrix &A) {
    assert(X.row() == A.row());
    assert(Y.row() == A.col());
    cblas_dger(CblasColMajor, X.row(), Y.row(), alpha, X.data(), 1, Y.data(), 1, A.data(), A.row());
}

/*
 * y = a Ax + b y
 */
void dgemv(scalar_t alpha, Matrix &A, Vector &x, scalar_t beta, Vector &y) {
    assert(A.col() == x.row());
    assert(A.row() == y.row());
    cblas_dgemv(CblasColMajor, CblasNoTrans, A.row(), A.col(), alpha, A.data(), A.row(), x.data(), 1, beta,
                y.data(), 1);
}

/*
 * y = a Atx + b y
 */
void dgemv_t(scalar_t alpha, Matrix &A, Vector &x, scalar_t beta, Vector &y) {
    assert(A.row() == x.row());
    assert(A.col() == y.row());
    cblas_dgemv(CblasColMajor, CblasTrans, A.row(), A.col(), alpha, A.data(), A.row(), x.data(), 1, beta,
                y.data(), 1);
}

/*
 * hadamard product
 */
void dsbmv(scalar_t alpha, Vector &A, Vector &x, scalar_t beta, Vector &y) {
    assert(A.row() == x.row());
    assert(y.row() == x.row());
    cblas_dsbmv(CblasColMajor, CblasLower, A.row(), 0, alpha, A.data(), 1, x.data(), 1, beta, y.data(), 1);
}

/*
 * unsafe case
 */
void dsbmv(scalar_t alpha, Vector &A, Vector &x, scalar_t beta, scalar_t *y) {
    assert(A.row() == x.row());
    cblas_dsbmv(CblasColMajor, CblasLower, A.row(), 0, alpha, A.data(), 1, x.data(), 1, beta, y, 1);
}

scalar_t nrm2(Vector &x) {
    assert(x.row() > 0);
    return cblas_dnrm2(x.row(), x.data(), 1);
}

scalar_t ddot(Vector &x, Vector &y) {
    assert(x.row() == y.row());
    return cblas_ddot(x.row(), x.data(), 1, y.data(), 1);
}
