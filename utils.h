//
// Created by lurker on 2/10/17.
//

#ifndef IBIM_OCTREE_UTILS_H
#define IBIM_OCTREE_UTILS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <utility>
#include <functional>
#include <chrono>
#include <iomanip>
#include "cblas.h"

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

#ifdef RUN_OMP
#include "omp.h"
#endif

#define EPS 1e-12
#define DIM 3 /* 3D FMM, do not change */

#define SQR(X) ((X)*(X))

#ifdef DISP
#define RUN(s, func){ \
std::chrono::steady_clock::time_point begin =std::chrono::steady_clock::now(); \
func;\
std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();\
std::cout << std::setw(15)<< s << " "  << std::setprecision(5) << std::setw(8) << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0 << " seconds"<<std::endl;\
}
#else
#define RUN(s, func){\
func;\
}
#endif

typedef int index_t;
typedef double scalar_t;
typedef bool bool_t;

using std::unordered_set;
using std::vector;
using std::queue;
using std::unordered_map;
using std::unordered_set;
using std::max;
using std::min;
using std::sqrt;
using std::acos;
using std::atan;

#endif //IBIM_OCTREE_UTILS_H
