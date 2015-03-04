# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([Grid], [1.0], [paboyle@ph.ed.ac.uk])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([Grid.h])
AC_CONFIG_HEADERS([Grid_config.h])

# Checks for programs.
AC_PROG_CXX
AC_OPENMP

# Checks for libraries.
AX_GCC_VAR_ATTRIBUTE(aligned)

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

AC_ARG_ENABLE([simd],[AC_HELP_STRING([--enable-simd=AVX|AVX2|SSE|SSE2|AVX512],[Select instructions])],[ac_SIMD=${enable_simd}],[ac_SIMD=AVX2])

case ${ac_SIMD} in
     SSE)
       echo Configuring for SSE1
       AC_DEFINE([SSE],[1],[SSE] )
     ;;
     SSE2)
       echo Configuring for SSE2
       AC_DEFINE([SSE2],[1],[SSE2] )
     ;;
     AVX)
       echo Configuring for AVX
       AC_DEFINE([AVX],[1],[AVX] )
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

AC_OUTPUT
