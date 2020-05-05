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

/* NVCC save and restore compile environment*/
#ifdef __NVCC__
#pragma push
#pragma diag_suppress code_is_unreachable
#pragma push_macro("__CUDA_ARCH__")
#pragma push_macro("__NVCC__")
#pragma push_macro("__CUDACC__")
#undef __NVCC__
#undef __CUDACC__
#undef __CUDA_ARCH__
#define __NVCC__REDEFINE__
#endif 

/* SYCL save and restore compile environment*/
#ifdef __SYCL_DEVICE_ONLY__  
#pragma push
#pragma push_macro("__SYCL_DEVICE_ONLY__")
#undef __SYCL_DEVICE_ONLY__
#undef EIGEN_USE_SYCL
#define EIGEN_DONT_VECTORIZE
#endif


#include <Grid/Eigen/Dense>
#include <Grid/Eigen/unsupported/CXX11/Tensor>

/* NVCC restore */
#ifdef __NVCC__REDEFINE__
#pragma pop_macro("__CUDACC__")
#pragma pop_macro("__NVCC__")
#pragma pop_macro("__CUDA_ARCH__")
#pragma pop
#endif

/*SYCL restore*/
#ifdef __SYCL__REDEFINE__
#pragma pop_macro("__SYCL_DEVICE_ONLY__")
#pragma pop
#endif

#if defined __GNUC__
#pragma GCC diagnostic pop
#endif


