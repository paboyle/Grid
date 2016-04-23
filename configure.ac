#                         -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
# Project Grid package  
# 
# Time-stamp: <2015-07-10 17:46:21 neo>

AC_PREREQ([2.63])
AC_INIT([Grid], [1.0], [paboyle@ph.ed.ac.uk])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE(subdir-objects)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([lib/Grid.h])
AC_CONFIG_HEADERS([lib/Config.h])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_MSG_NOTICE([

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Configuring $PACKAGE v$VERSION  for $host
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
])

# Checks for programs.
AC_LANG(C++)
AC_PROG_CXX
AC_OPENMP
AC_PROG_RANLIB
#AX_CXX_COMPILE_STDCXX_11(noext, mandatory)
AX_EXT

# Checks for libraries.
#AX_GCC_VAR_ATTRIBUTE(aligned)

# Checks for header files.
AC_CHECK_HEADERS(stdint.h)
AC_CHECK_HEADERS(mm_malloc.h)
AC_CHECK_HEADERS(malloc/malloc.h)
AC_CHECK_HEADERS(malloc.h)
AC_CHECK_HEADERS(endian.h)
AC_CHECK_HEADERS(execinfo.h)
AC_CHECK_HEADERS(gmp.h)
AC_CHECK_DECLS([ntohll],[], [], [[#include <arpa/inet.h>]])
AC_CHECK_DECLS([be64toh],[], [], [[#include <arpa/inet.h>]])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
echo
echo Checking libraries 
echo :::::::::::::::::::::::::::::::::::::::::::

AC_CHECK_FUNCS([gettimeofday])


#
# SIMD instructions selection
#

AC_ARG_ENABLE([simd],[AC_HELP_STRING([--enable-simd=SSE4|AVX|AVXFMA4|AVX2|AVX512|IMCI],\
	[Select instructions to be SSE4.0, AVX 1.0, AVX 2.0+FMA, AVX 512, IMCI])],\
	[ac_SIMD=${enable_simd}],[ac_SIMD=DEBUG])

supported=no

ac_ZMM=no;

case ${ac_SIMD} in
     SSE4)
       echo Configuring for SSE4
       AC_DEFINE([SSE4],[1],[SSE4 Intrinsics] )
       if test x"$ax_cv_support_ssse3_ext" = x"yes"; then  dnl minimal support for SSE4
         supported=yes
       else
  	AC_MSG_WARN([Your processor does not support SSE4 instructions])
       fi
     ;;
     AVX)
       echo Configuring for AVX
       AC_DEFINE([AVX1],[1],[AVX Intrinsics] )
       if test x"$ax_cv_support_avx_ext" = x"yes"; then  dnl minimal support for AVX
       supported=yes			  
       else
       	AC_MSG_WARN([Your processor does not support AVX instructions])
       fi
     ;;
     AVXFMA4)
       echo Configuring for AVX
       AC_DEFINE([AVXFMA4],[1],[AVX Intrinsics with FMA4] )
       if test x"$ax_cv_support_avx_ext" = x"yes"; then  dnl minimal support for AVX
       supported=yes			  
       else
       	AC_MSG_WARN([Your processor does not support AVX instructions])
       fi
     ;;
     AVX2)
       echo Configuring for AVX2
       AC_DEFINE([AVX2],[1],[AVX2 Intrinsics] )
       if test x"$ax_cv_support_avx2_ext" = x"yes"; then  dnl minimal support for AVX2
       supported=yes
       else
       AC_MSG_WARN([Your processor does not support AVX2 instructions])
       fi
     ;;
     AVX512)
       echo Configuring for AVX512 
       AC_DEFINE([AVX512],[1],[AVX512 Intrinsics for Knights Landing] )
       supported="cross compilation"
       ac_ZMM=yes;
     ;;
     IMCI)
       echo Configuring for IMCI
       AC_DEFINE([IMCI],[1],[IMCI Intrinsics for Knights Corner] )
       supported="cross compilation"
       ac_ZMM=no;
     ;;
     NEONv8)
       echo Configuring for experimental ARMv8a support 
       AC_DEFINE([NEONv8],[1],[NEON ARMv8 Experimental support ] )
       supported="cross compilation"
     ;;
     DEBUG)
       echo Configuring without SIMD support - only for compiler DEBUGGING!
       AC_DEFINE([EMPTY_SIMD],[1],[EMPTY_SIMD only for DEBUGGING] )
      ;;     
     *)
     AC_MSG_ERROR([${ac_SIMD} flag unsupported as --enable-simd option\nRun ./configure --help for the list of options]); 
     ;;
esac

case ${ac_ZMM} in
yes)
	echo Enabling ZMM source code
;;
no)
	echo Disabling ZMM source code
;;
esac

AM_CONDITIONAL(BUILD_ZMM,[ test "X${ac_ZMM}X" == "XyesX" ])

AC_ARG_ENABLE([precision],[AC_HELP_STRING([--enable-precision=single|double],[Select default word size of Real])],[ac_PRECISION=${enable_precision}],[ac_PRECISION=double])
case ${ac_PRECISION} in
     single)
       echo default precision is single
       AC_DEFINE([GRID_DEFAULT_PRECISION_SINGLE],[1],[GRID_DEFAULT_PRECISION is SINGLE] )
     ;;
     double)
       echo default precision is double
       AC_DEFINE([GRID_DEFAULT_PRECISION_DOUBLE],[1],[GRID_DEFAULT_PRECISION is DOUBLE] )
     ;;
esac

#
# Comms selection
#

AC_ARG_ENABLE([comms],[AC_HELP_STRING([--enable-comms=none|mpi],[Select communications])],[ac_COMMS=${enable_comms}],[ac_COMMS=none])

case ${ac_COMMS} in
     none)
       echo Configuring for NO communications
       AC_DEFINE([GRID_COMMS_NONE],[1],[GRID_COMMS_NONE] )
     ;;
     mpi)
       echo Configuring for MPI communications
       AC_DEFINE([GRID_COMMS_MPI],[1],[GRID_COMMS_MPI] )
     ;;
     shmem)
       echo Configuring for SHMEM communications
       AC_DEFINE([GRID_COMMS_SHMEM],[1],[GRID_COMMS_SHMEM] )
     ;;
     *)
     AC_MSG_ERROR([${ac_COMMS} unsupported --enable-comms option]); 
     ;;
esac

AM_CONDITIONAL(BUILD_COMMS_SHMEM,[ test "X${ac_COMMS}X" == "XshmemX" ])
AM_CONDITIONAL(BUILD_COMMS_MPI,[ test "X${ac_COMMS}X" == "XmpiX" ])
AM_CONDITIONAL(BUILD_COMMS_NONE,[ test "X${ac_COMMS}X" == "XnoneX" ])

#
# RNG selection
#
AC_ARG_ENABLE([rng],[AC_HELP_STRING([--enable-rng=ranlux48|mt19937],\
	[Select Random Number Generator to be used])],\
	[ac_RNG=${enable_rng}],[ac_RNG=ranlux48])
case ${ac_RNG} in
     ranlux48)
     AC_DEFINE([RNG_RANLUX],[1],[RNG_RANLUX] )
     ;;
     mt19937)
     AC_DEFINE([RNG_MT19937],[1],[RNG_MT19937] )
     ;;
     *)
     AC_MSG_ERROR([${ac_RNG} unsupported --enable-rng option]); 
     ;;
esac
#
# Chroma regression tests
#
AC_ARG_ENABLE([chroma],[AC_HELP_STRING([--enable-chroma],[Expect chroma compiled under c++11 ])],ac_CHROMA=yes,ac_CHROMA=no)

case ${ac_CHROMA} in
     yes)
       echo Enabling tests regressing to Chroma
     ;;
     no)
       echo Disabling tests regressing to Chroma
     ;;
     *)
     AC_MSG_ERROR([${ac_CHROMA} unsupported --enable-chroma option]); 
     ;;
esac

AM_CONDITIONAL(BUILD_CHROMA_REGRESSION,[ test "X${ac_CHROMA}X" == "XyesX" ])

#
# Lapack
#
AC_ARG_ENABLE([lapack],[AC_HELP_STRING([--enable-lapack],[Enable lapack yes/no ])],[ac_LAPACK=${enable_lapack}],[ac_LAPACK=no])

case ${ac_LAPACK} in
     yes)
       echo Enabling lapack
     ;;
     no)
       echo Disabling lapack
     ;;
     *)
       echo Enabling lapack at ${ac_LAPACK}
     ;;
esac

AM_CONDITIONAL(USE_LAPACK,[ test "X${ac_LAPACK}X" != "XnoX" ])
AM_CONDITIONAL(USE_LAPACK_LIB,[ test "X${ac_LAPACK}X" != "XyesX" ])

###################################################################
# Checks for doxygen support
# if present enables the "make doxyfile" command
#echo
#echo Checking doxygen support 
#echo :::::::::::::::::::::::::::::::::::::::::::
#AC_PROG_DOXYGEN

#if test -n "$DOXYGEN"
#then
#AC_CONFIG_FILES([docs/doxy.cfg])
#fi

echo
echo Creating configuration files
echo :::::::::::::::::::::::::::::::::::::::::::
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(tests/Makefile)
AC_CONFIG_FILES(tests/qdpxx/Makefile)
AC_CONFIG_FILES(benchmarks/Makefile)
AC_OUTPUT


echo "
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary of configuration for $PACKAGE v$VERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following features are enabled:

- architecture (build)          : $build_cpu
- os (build)                    : $build_os
- architecture (target)         : $target_cpu
- os (target)                   : $target_os
- build DOXYGEN documentation   : `if test "x$enable_doc" = xyes; then echo yes; else echo no; fi`
- graphs and diagrams           : `if test "x$enable_dot" = xyes; then echo yes; else echo no; fi`
- Supported SIMD flags          : $SIMD_FLAGS
----------------------------------------------------------
- enabled simd support          : ${ac_SIMD}   (config macro says supported: $supported )
- communications type           : ${ac_COMMS}
- default precision             : ${ac_PRECISION}
- RNG choice                    : ${ac_RNG} 
- LAPACK	                : ${ac_LAPACK} 


"
