#include <Grid/GridCore.h>
#pragma once
// Force Eigen to use MKL if Grid has been configured with --enable-mkl
#ifdef USE_MKL
#define EIGEN_USE_MKL_ALL
#endif

#if defined __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifdef __NVCC__

#pragma push
#pragma diag_suppress code_is_unreachable

#define __NVCC__REDEFINE__

#undef __NVCC__
#undef __CUDACC__

#ifdef __CUDA_ARCH__

#define __CUDA_ARCH__REDEFINE__ 1
#define __CUDA_ARCH_SAVE__ __CUDA_ARCH__

#undef __CUDA_ARCH__

#endif

#endif

#include <Grid/Eigen/Dense>

#ifdef __NVCC__REDEFINE__

#pragma pop

#define __NVCC__
#define __CUDACC__

#ifdef __CUDA_ARCH__REDEFINE__

#define __CUDA_ARCH__ __CUDA_ARCH_SAVE__

#endif
#endif

#if defined __GNUC__
#pragma GCC diagnostic pop
#endif
