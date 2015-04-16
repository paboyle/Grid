# Process this file with autoconf to produce a configure script.
AC_INIT([Grid], [1.0], [paboyle@ph.ed.ac.uk])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([Grid.h])
AC_CONFIG_HEADERS([Grid_config.h])

# Checks for programs.
AC_PROG_CXX
AC_OPENMP
AC_PROG_RANLIB

# Checks for libraries.
#AX_GCC_VAR_ATTRIBUTE(aligned)

# Checks for header files.
AC_CHECK_HEADERS(stdint.h)
AC_CHECK_HEADERS(malloc/malloc.h)
AC_CHECK_HEADERS(malloc.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday])

AC_ARG_ENABLE([simd],[AC_HELP_STRING([--enable-simd=SSE|AVX|AVX2|AVX512],[Select instructions])],[ac_SIMD=${enable_simd}],[ac_SIMD=AVX2])

case ${ac_SIMD} in
     SSE4)
       echo Configuring for SSE4
       AC_DEFINE([SSE4],[1],[SSE4] )
     ;;
     AVX)
       echo Configuring for AVX
       AC_DEFINE([AVX1],[1],[AVX] )
     ;;
     AVX2)
       echo Configuring for AVX2
       AC_DEFINE([AVX2],[1],[AVX2] )
     ;;
     AVX512)
       echo Configuring for AVX512
       AC_DEFINE([AVX512],[1],[AVX512] )
     ;;
     *)
     AC_MSG_ERROR([${ac_SIMD} unsupported --enable-simd option]); 
     ;;
esac


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
     *)
     AC_MSG_ERROR([${ac_COMMS} unsupported --enable-comms option]); 
     ;;
esac

AM_CONDITIONAL(BUILD_COMMS_MPI,[ test "X${ac_COMMS}X" == "XmpiX" ])
AM_CONDITIONAL(BUILD_COMMS_NONE,[ test "X${ac_COMMS}X" == "XnoneX" ])


AC_CONFIG_FILES(Makefile)
AC_OUTPUT
