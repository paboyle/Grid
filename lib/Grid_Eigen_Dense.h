#include <Grid/GridCore.h>
#pragma once
#if defined __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#ifdef GRID_NVCC
#pragma push
#pragma diag_suppress code_is_unreachable
#undef __NVCC__
#undef __CUDACC__
#ifdef __CUDA_ARCH__
#define __CUDA_ARCH__REDEFINE__
#undef __CUDA_ARCH__
#endif
#endif

#include <Grid/Eigen/Dense>

#ifdef GRID_NVCC
#pragma pop
#define __NVCC__
#define __CUDACC__

#ifdef __CUDA_ARCH__REDEFINE__
#define __CUDA_ARCH__
#endif
#endif

#if defined __GNUC__
#pragma GCC diagnostic pop
#endif
