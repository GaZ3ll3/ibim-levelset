//
// Created by lurker on 3/22/17.
//

#ifndef LEVELSET_BBFMM_H
#define LEVELSET_BBFMM_H


#if !defined __extern_always_inline && defined __clang__
# if defined __GNUC_STDC_INLINE__ || defined __GNUC_GNU_INLINE__
#  define __extern_inline extern __inline __attribute__ ((__gnu_inline__))
#  define __extern_always_inline \
  extern __always_inline __attribute__ ((__gnu_inline__))
# else
#  define __extern_inline extern __inline
#  define __extern_always_inline extern __always_inline
# endif
#endif


#include "utils.h"
#include "blas_wraper.h"

namespace bbfmm {
    class basePoint {
    public:
        scalar_t x;
        scalar_t y;
        scalar_t z;

        basePoint() : x(0.), y(0.), z(0.) {}

        basePoint(scalar_t _x, scalar_t _y, scalar_t _z) : x(_x), y(_y), z(_z) {}

        virtual ~basePoint() {}

        bool operator>=(const basePoint &a) {
            return (x >= a.x - EPS) && (y >= a.y - EPS) && (z >= a.z - EPS);
        }

        bool operator<=(const basePoint &a) {
            return (x <= a.x + EPS) && (y <= a.y + EPS) && (z <= a.z + EPS);
        }

        bool operator==(const basePoint &a) {
            return fabs(x - a.x) < EPS && fabs(y - a.y) < EPS && fabs(z - a.z) < EPS;
        }
    };


    class point : public basePoint {
    public:
        /*
         * additional variables
         */
        int triangleId;
        int sphereId;

        point() : basePoint() {
            triangleId = -1;
            sphereId = -1;
        }

        point(scalar_t _x, scalar_t _y, scalar_t _z) : basePoint(_x, _y, _z) {
            triangleId = -1;
            sphereId = -1;
        }

        ~point() {}
    };


    class baseNode {
    public:
        index_t parent;
        index_t child[8];

        unordered_set <index_t> uList;
        unordered_set <index_t> vList;
        unordered_set <index_t> wList;
        unordered_set <index_t> xList;

        index_t nLevel;

        index_t nUList;
        index_t nVList;
        index_t nWList;
        index_t nXList;

        index_t nodeIndex;
        point center;
        point radius;

        index_t nSource;
        index_t nTarget;

        vector <index_t> sourceIndex;
        vector <index_t> targetIndex;

        bool isLeaf;
        bool isEmpty;
        bool chargeComputed;

        baseNode(index_t level, index_t index) {
            parent = -1;
            for (index_t i = 0; i < 8; ++i) {
                child[i] = -1;
            }
            nLevel = level;
            nodeIndex = index;
            isLeaf = false;
            isEmpty = false;
            chargeComputed = false;
            nSource = 0;
            nTarget = 0;
            nUList = 0;
            nVList = 0;
            nWList = 0;
            nXList = 0;

        }

        virtual ~baseNode() {}

    };


    class node : public baseNode {
    public:
        node(index_t level, index_t index) : baseNode(level, index) {

        }

        ~node() {}

        /*
         *  some internal members
         */
        vector <point> scaledCnode;
        Vector potential;
        Vector nodeCharge;
        Vector charge;
        Vector nodePotential;
        Matrix R;
        Matrix L;
    };


    class tree {
    public:
        vector <node> dict;
        index_t maxId;
        index_t root;
        index_t nSource;
        index_t nTarget;

        index_t rank;
        index_t maxLevel;

        vector <point> sourceTree;
        vector <point> targetTree;

        point center;
        point radius;

        tree() {
            maxId = -1;
            root = -1;
            nSource = 0;
            nTarget = 0;
            rank = 0;
            maxLevel = 0;
        }

        ~tree() {

        }

        void
        populate(vector <point> &_source, vector <point> &_target, index_t _nSource, index_t _nTarget, index_t _rank,
                 index_t _maxLevel);

        void output(std::string file);

    protected:
        void getCenterRadius(vector <point> &_source);

        void assignChildren(index_t _id, index_t _maxLevel);

        void buildTree();

        void buildNode(index_t _id, point &min_p, point &max_p);

        index_t findNode(index_t _id, point &p);

        bool isAdjacent(index_t _aId, index_t _bId);

    };

    class kernel {
    public:
        tree t;
        Vector chargeTree;
        std::function<scalar_t(point &, point &)> eval;
        index_t rank;

        Matrix R[8];

        index_t nChebyshev;
        Vector chebyNode;
        Matrix tNode;

        kernel() {
            nChebyshev = 0;
            rank = 0;
        }

        ~kernel() {}

        void
        initialize(index_t _nChebyshev, vector<point> &_source, vector<point> &_target, Vector &_charge,
                   index_t _nSource,
                   index_t _nTarget, index_t _rank, index_t _maxLevel) {
            // populate the kd-tree.
            t.populate(_source, _target, _nSource, _nTarget, _rank, _maxLevel);

            nChebyshev = _nChebyshev;

            chargeTree = _charge;

            // nChebyshev^3 nodes are used for interpolation.
            rank = nChebyshev * nChebyshev * nChebyshev;

            chebyNode = Vector(nChebyshev);
            getStandardChebyNodes(nChebyshev, chebyNode);

            tNode = Matrix(nChebyshev, nChebyshev);
            getStandardChebyPoly(nChebyshev, nChebyshev, chebyNode, tNode);

            getTransfer(nChebyshev, chebyNode, tNode, R);

        }

        void run(Vector &potentialMatrix) {
#ifdef RUN_OMP
#pragma omp parallel
#endif
            {
#ifdef RUN_OMP
#pragma omp single
#endif
                RUN("up-pass", upPass(0));
            }
#ifdef RUN_OMP
#pragma omp taskwait
#endif
            potentialMatrix = Vector(t.nTarget);
#ifdef RUN_OMP
#pragma omp parallel
#endif
            {
#ifdef RUN_OMP
#pragma  omp single
#endif
                RUN("down-pass", downPass(0, potentialMatrix));
            }
#ifdef RUN_OMP
#pragma omp taskwait
#endif
        }


        void getStandardChebyNodes(index_t _nChebyshev, Vector &_chebyNode) {
            assert(_chebyNode.row() == nChebyshev);
            for (index_t i = 0; i < _nChebyshev; ++i) {
                _chebyNode(i) = -cos((i + 0.5) * M_PI / _nChebyshev);
            }
        }

        void getStandardChebyPoly(index_t _nChebyPoly, index_t _N, Vector &_x, Matrix &_T) {
            assert(_T.row() == _N);
            assert(_T.col() == _nChebyPoly);

            setValue(_T, 0);

            Vector ones(_N);
            setValue(ones, 1.0);
            _T.setColumn(0, ones);

            if (_nChebyPoly > 1) {
                _T.setColumn(1, _x);
                for (index_t i = 2; i < _nChebyPoly; ++i) {
                    /*
                     * only copy pointers
                     */
                    Vector T1(_N, false, _T.column(i - 1));
                    Vector T2(_N, false, _T.column(i - 2));

                    dsbmv(2.0, _x, T1, 0., _T.column(i));
                    daxpy(-1.0, T2, _T.column(i));
                }
            }
        }

        void getTransferFromParentChebyshevToChildrenChebyshev(index_t _nChebyshev, Vector &_chebyNode, Matrix &_tNode,
                                                               Matrix &_transfer) {
            Vector childChebyNode(2 * _nChebyshev);
            Vector T1(_nChebyshev);
            setValue(T1, -0.5);
            daxpy(0.5, _chebyNode, T1);
            memcpy(childChebyNode.data(), T1.data(), _nChebyshev * sizeof(scalar_t));

            Vector T2(_nChebyshev);
            setValue(T2, 0.5);
            daxpy(0.5, _chebyNode, T2);
            memcpy(childChebyNode.data() + _nChebyshev, T2.data(), _nChebyshev * sizeof(scalar_t));

            getStandardChebyPoly(_nChebyshev, 2 * _nChebyshev, childChebyNode, _transfer);

            Matrix T3(2 * _nChebyshev, _nChebyshev);
            setValue(T3, -1.0);

            dgemm_t(2.0, _transfer, _tNode, 1.0, T3);
            dscal(1.0 / _nChebyshev, T3);
            _transfer = T3;
        }

        void getTransfer(index_t _nChebyshev, Vector &_chebyNode, Matrix &_tNode, Matrix *R) {
            Matrix S(2 * _nChebyshev, _nChebyshev);
            getTransferFromParentChebyshevToChildrenChebyshev(_nChebyshev, _chebyNode, _tNode, S);

            Matrix Transfer[2];
            Transfer[0].resize(_nChebyshev, _nChebyshev);
            Transfer[1].resize(_nChebyshev, _nChebyshev);

            setBlock(Transfer[0], S, 0, 0, _nChebyshev, _nChebyshev);
            setBlock(Transfer[1], S, _nChebyshev, 0, _nChebyshev, _nChebyshev);

            index_t _rank = _nChebyshev * _nChebyshev * _nChebyshev;
            for (index_t i = 0; i < 8; ++i) {
                R[i].resize(_rank, _rank);
            }

            // follow bit representaion.
            for (index_t i = 0; i < _nChebyshev; ++i) {
                for (index_t j = 0; j < _nChebyshev; ++j) {
                    for (index_t k = 0; k < _nChebyshev; ++k) {
                        for (index_t l = 0; l < _nChebyshev; ++l) {
                            for (index_t m = 0; m < _nChebyshev; ++m) {
                                for (index_t n = 0; n < _nChebyshev; ++n) {
                                    for (index_t id = 0; id < 8; ++id) {
                                        index_t bit[3];
                                        bit[0] = (id >> 0) & 1;
                                        bit[1] = (id >> 1) & 1;
                                        bit[2] = (id >> 2) & 1;
                                        R[id](i * _nChebyshev * _nChebyshev + j * _nChebyshev + k,
                                              l * _nChebyshev * _nChebyshev + m * _nChebyshev + n) =
                                                Transfer[bit[2]](i, l) * Transfer[bit[1]](j, m) *
                                                Transfer[bit[0]](k, n);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        void getScaledChebyNode(index_t _nChebyNode, Vector &_chebyNode, point &center, point &radius,
                                vector<point> &_scaledCnode) {
            for (index_t i = 0; i < _nChebyNode; ++i) {
                _scaledCnode.push_back(point(center.x + radius.x * _chebyNode(i),
                                             center.y + radius.y * _chebyNode(i),
                                             center.z + radius.z * _chebyNode(i)));
            }
        }

        void getCharge(index_t rootId) {
            node &n = t.dict[rootId];
            if (n.chargeComputed) {
                return;
            } else {
                n.chargeComputed = true;
                n.charge.resize(n.nSource);
                for (index_t k = 0; k < n.nSource; ++k) {
                    n.charge(k) = chargeTree(n.sourceIndex[k]);
                }
            }
        }

        void
        getTransferParentToChildren(index_t _nChebyNode, vector<point> &_tree, vector<index_t> &_index, point &_center,
                                    point &_radius,
                                    Vector &_chebyNode, Matrix &_tNode, Matrix &R) {

            index_t N = (index_t) _index.size();
            Vector standlocation[3];
            standlocation[0].resize(N);
            standlocation[1].resize(N);
            standlocation[2].resize(N);
            for (index_t i = 0; i < N; ++i) {
                standlocation[0](i) = (_tree[_index[i]].x - _center.x) / _radius.x;
                standlocation[1](i) = (_tree[_index[i]].y - _center.y) / _radius.y;
                standlocation[2](i) = (_tree[_index[i]].z - _center.z) / _radius.z;
            }

            Matrix Transfer[3];
            for (index_t k = 0; k < 3; ++k) {
                Transfer[k].resize(N, _nChebyNode);
                getStandardChebyPoly(_nChebyNode, N, standlocation[k], Transfer[k]);

                Matrix T3(N, _nChebyNode);
                setValue(T3, -1.0);
                dgemm_t(2.0, Transfer[k], _tNode, 1.0, T3);
                dscal(1.0 / _nChebyNode, T3);
                Transfer[k] = T3;

            }

            index_t _rank = _nChebyNode * _nChebyNode * _nChebyNode;
            R.resize(N, _rank);
            for (index_t k = 0; k < N; ++k) {
                for (index_t i = 0; i < _nChebyNode; ++i) {
                    for (index_t j = 0; j < _nChebyNode; ++j) {
                        for (index_t l = 0; l < _nChebyNode; ++l) {
                            R(k, l * _nChebyNode * _nChebyNode + j * _nChebyNode + i) =
                                    Transfer[0](k, i) * Transfer[1](k, j) * Transfer[2](k, l);
                        }
                    }
                }
            }
        }

        void kernelEval(vector<point> &_source, vector<point> &_target, Matrix &K) {
            K.resize((index_t) _target.size(), (index_t) _source.size());
            for (index_t _s = 0; _s < _source.size(); ++_s) {
                for (index_t _t = 0; _t < _target.size(); ++_t) {
                    K(_t, _s) = this->eval(_source[_s], _target[_t]);
                }
            }
        }

        void kernelEvalIndex(vector<index_t> &_sourceIndex, vector<index_t> &_targetIndex, Matrix &K) {
            K.resize((index_t) _targetIndex.size(), (index_t) _sourceIndex.size());
            for (index_t _s = 0; _s < _sourceIndex.size(); ++_s) {
                for (index_t _t = 0; _t < _targetIndex.size(); ++_t) {
                    K(_t, _s) = this->eval(
                            this->t.sourceTree[_sourceIndex[_s]],
                            this->t.targetTree[_targetIndex[_t]]
                    );
                }
            }
        }

        void kernelEvalChebyshev(index_t _M, vector<point> &_xv, index_t _N, vector<point> &_yv, Matrix &K) {
            vector<point> sourceVec;
            vector<point> targetVec;
            K.resize(_M * _M * _M, _N * _N * _N);
            for (index_t k = 0; k < _M; k++) {
                for (index_t j = 0; j < _M; j++) {
                    for (index_t i = 0; i < _M; i++) {
                        point np(_xv[i].x, _xv[j].y, _xv[k].z);
                        sourceVec.push_back(np);
                    }
                }
            }

            for (index_t k = 0; k < _N; k++) {
                for (index_t j = 0; j < _N; j++) {
                    for (index_t i = 0; i < _N; i++) {
                        point np(_yv[i].x, _yv[j].y, _yv[k].z);
                        targetVec.push_back(np);
                    }
                }
            }

            kernelEval(sourceVec, targetVec, K);
        }

        void reset(index_t rootId = 0) {
            if (rootId < 0)
                return;
            node &n = t.dict[rootId];

            n.chargeComputed = false;
            n.scaledCnode.clear();
            setValue(n.nodeCharge, 0.);
            setValue(n.nodePotential, 0.);
            setValue(n.potential, 0.);
            setValue(n.R, 0.);
            setValue(n.L, 0.);

            for (index_t i = 0; i < 8; ++i) {
                reset(n.child[i]);
            }
        }


        void upPass(index_t rootId = 0) {
            node &n = t.dict[rootId];
            n.scaledCnode.clear();
            n.nodeCharge.resize(rank);
            n.nodePotential.resize(rank);
            getScaledChebyNode(nChebyshev, chebyNode, n.center, n.radius, n.scaledCnode);

            if (n.isLeaf) {
                // lazy
                getCharge(rootId);
                getTransferParentToChildren(nChebyshev, t.sourceTree, n.sourceIndex, n.center, n.radius, chebyNode,
                                            tNode,
                                            n.R);
                getTransferParentToChildren(nChebyshev, t.targetTree, n.targetIndex, n.center, n.radius, chebyNode,
                                            tNode,
                                            n.L);
                /*
                 * in case leaf node has no points, which causes DGEMV error.
                 */
                if (n.R.row() != 0) dgemv_t(1.0, n.R, n.charge, 1.0, n.nodeCharge);
            } else {
                for (index_t i = 0; i < 8; ++i) {
#ifdef RUN_OMP
#pragma omp task shared(n) firstprivate(i)
#endif
                    upPass(n.child[i]);
                }
#ifdef RUN_OMP
#pragma omp taskwait
#endif
                for (index_t i = 0; i < 8; ++i) {
                    if (!t.dict[n.child[i]].isEmpty) {
                        dgemv_t(1.0, R[i], t.dict[n.child[i]].nodeCharge, 1.0, n.nodeCharge);
                    }
                }
            }
        }

        void downPass(index_t rootId, Vector &potential) {
            node &n = t.dict[rootId];
            Matrix K;

            Vector temp;
            if (n.parent != -1) {
                /*
                 * V list
                 */
                for (index_t i : n.vList) {
                    if (!t.dict[i].isEmpty) {
                        kernelEvalChebyshev(nChebyshev, t.dict[i].scaledCnode, nChebyshev, n.scaledCnode, K);
                        dgemv(1.0, K, t.dict[i].nodeCharge, 1.0, n.nodePotential);
                    }
                }
                /*
                 * X List
                 */
                for (index_t i : n.xList) {
                    if (!t.dict[i].isEmpty) {
                        kernelEvalChebyshev(nChebyshev, t.dict[i].scaledCnode, nChebyshev, n.scaledCnode, K);
                        dgemv(1.0, K, t.dict[i].nodeCharge, 1.0, n.nodePotential);
                    }
                }

                /*
                 * L2L
                 */
                node &p = t.dict[n.parent];
                dgemv(1.0, this->R[n.nodeIndex], p.nodePotential, 1.0, n.nodePotential);

            }

            if (n.isLeaf && n.nTarget != 0) {
                n.potential.resize(n.nTarget);

                /*
                 * U List
                 */
                for (index_t i : n.uList) {
                    if (!t.dict[i].isEmpty) {
                        getCharge(i);
                        kernelEvalIndex(t.dict[i].sourceIndex, n.targetIndex, K);
                        dgemv(1.0, K, t.dict[i].charge, 1.0, n.potential);
                    }
                }

                /*
                 * W List
                 */

                for (index_t i : n.wList) {
                    if (!t.dict[i].isEmpty) {
                        getCharge(i);
                        kernelEvalIndex(t.dict[i].sourceIndex, n.targetIndex, K);
                        dgemv(1.0, K, t.dict[i].charge, 1.0, n.potential);
                    }
                }

                /*
                 * L2T
                 */
                dgemv(1.0, n.L, n.nodePotential, 1.0, n.potential);

                /*
                 * Finalize, caution:
                 *
                 * omp should be fine here, because no two threads will write to the same place at the same time.
                 */
                for (index_t i = 0; i < n.nTarget; i++) {
                    potential(n.targetIndex[i]) += n.potential(i);
                }
            }

            if (!n.isLeaf) {
                for (index_t i = 0; i < 8; ++i) {
#ifdef RUN_OMP
#pragma omp task shared(n, potential) firstprivate(i)
#endif
                    downPass(n.child[i], potential);
                }
#ifdef RUN_OMP
#pragma omp taskwait
#endif
            }
        }
    };
}

#endif //LEVELSET_BBFMM_H
