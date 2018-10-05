#pragma once
// Force Eigen to use MKL if Grid has been configured with --enable-mkl
#ifdef USE_MKL
#define EIGEN_USE_MKL_ALL
#endif

#if defined __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#include <Grid/Eigen/Dense>
#if defined __GNUC__
#pragma GCC diagnostic pop
#endif
