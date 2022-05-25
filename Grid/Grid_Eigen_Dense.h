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
#if (__CUDACC_VER_MAJOR__ >= 11) && (__CUDACC_VER_MINOR__ >= 5)
#pragma nv_diag_suppress code_is_unreachable
#else
#pragma diag_suppress code_is_unreachable
#endif
#pragma push_macro("__CUDA_ARCH__")
#pragma push_macro("__NVCC__")
#pragma push_macro("__CUDACC__")
#undef __CUDA_ARCH__
#undef __NVCC__
#undef __CUDACC__
#define __NVCC__REDEFINE__
#endif 

/* SYCL save and restore compile environment*/
#ifdef GRID_SYCL
#pragma push
#pragma push_macro("__SYCL_DEVICE_ONLY__")
#undef __SYCL_DEVICE_ONLY__
#define EIGEN_DONT_VECTORIZE
//#undef EIGEN_USE_SYCL
#define __SYCL__REDEFINE__
#endif

/* HIP save and restore compile environment*/
#ifdef GRID_HIP
#pragma push
#pragma push_macro("__HIP_DEVICE_COMPILE__")
#endif
#define EIGEN_NO_HIP

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

/*HIP restore*/
#ifdef __HIP__REDEFINE__
#pragma pop_macro("__HIP_DEVICE_COMPILE__")
#pragma pop
#endif

#if defined __GNUC__
#pragma GCC diagnostic pop
#endif


