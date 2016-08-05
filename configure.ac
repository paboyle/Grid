AC_PREREQ([2.63])
AC_INIT([Grid], [0.5.1-dev], [https://github.com/paboyle/Grid], [Grid])
AM_INIT_AUTOMAKE(subdir-objects)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([lib/Grid.h])
AC_CONFIG_HEADERS([lib/Config.h])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

############### Checks for programs
AC_LANG(C++)
: ${CXXFLAGS="-O3"}
AC_PROG_CXX
AC_OPENMP
AM_CXXFLAGS="$OPENMP_CXXFLAGS $AM_CXXFLAGS"
LT_INIT([disable-shared])

############### Checks for header files
AC_CHECK_HEADERS(stdint.h)
AC_CHECK_HEADERS(mm_malloc.h)
AC_CHECK_HEADERS(malloc/malloc.h)
AC_CHECK_HEADERS(malloc.h)
AC_CHECK_HEADERS(endian.h)
AC_CHECK_HEADERS(execinfo.h)
AC_CHECK_DECLS([ntohll],[], [], [[#include <arpa/inet.h>]])
AC_CHECK_DECLS([be64toh],[], [], [[#include <arpa/inet.h>]])

############### Checks for typedefs, structures, and compiler characteristics
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

############### Options
AC_ARG_WITH([gmp],
    [AS_HELP_STRING([--with-gmp=prefix],
    [try this for a non-standard install prefix of the GMP library])],
    [AM_CXXFLAGS="-I$with_gmp/include $AM_CXXFLAGS"]
    [AM_LDFLAGS="-L$with_gmp/lib" $AM_LDFLAGS])
AC_ARG_WITH([mpfr],
    [AS_HELP_STRING([--with-mpfr=prefix],
    [try this for a non-standard install prefix of the MPFR library])],
    [AM_CXXFLAGS="-I$with_mpfr/include $AM_CXXFLAGS"]
    [AM_LDFLAGS="-L$with_mpfr/lib $AM_LDFLAGS"])
AC_ARG_ENABLE([lapack],
    [AC_HELP_STRING([--enable-lapack=yes|no|prefix], [enable LAPACK])], 
    [ac_LAPACK=${enable_lapack}],[ac_LAPACK=no])
case ${ac_LAPACK} in
    no)
        ;;
    yes)
        AC_DEFINE([USE_LAPACK],[1],[use LAPACK]);;
    *)
        AM_CXXFLAGS="-I$ac_LAPACK/include $AM_CXXFLAGS"
        AM_LDFLAGS="-L$ac_LAPACK/lib $AM_LDFLAGS"
        AC_DEFINE([USE_LAPACK],[1],[use LAPACK])
esac

################ Get compiler informations
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
AX_COMPILER_VENDOR
AC_DEFINE_UNQUOTED([CXX_COMP_VENDOR],["$ax_cv_cxx_compiler_vendor"],
      [vendor of C++ compiler that will compile the code])
AX_GXX_VERSION
AC_DEFINE_UNQUOTED([GXX_VERSION],["$GXX_VERSION"],
      [version of g++ that will compile the code])

############### Checks for library functions
CXXFLAGS_CPY=$CXXFLAGS
LDFLAGS_CPY=$LDFLAGS
LIBS_CPY=$LIBS
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"
AC_CHECK_FUNCS([gettimeofday])
AC_CHECK_LIB([gmp],[__gmpf_init],
             [AC_CHECK_LIB([mpfr],[mpfr_init],
                 [AC_DEFINE([HAVE_LIBMPFR], [1], [Define to 1 if you have the `MPFR' library (-lmpfr).])]
                 [have_mpfr=true]
                 [LIBS="$LIBS -lmpfr"],
                 [AC_MSG_ERROR([MPFR library not found])])]
             [AC_DEFINE([HAVE_LIBGMP], [1], [Define to 1 if you have the `GMP' library (-lgmp).])]
             [have_gmp=true]
             [LIBS="$LIBS -lgmp"],
             [AC_MSG_WARN([**** GMP library not found, Grid can still compile but RHMC will not work ****])])

if test "${ac_LAPACK}x" != "nox"; then
    AC_CHECK_LIB([lapack],[LAPACKE_sbdsdc],[],
                 [AC_MSG_ERROR("LAPACK enabled but library not found")])
fi
CXXFLAGS=$CXXFLAGS_CPY
LDFLAGS=$LDFLAGS_CPY

############### SIMD instruction selection
AC_ARG_ENABLE([simd],[AC_HELP_STRING([--enable-simd=SSE4|AVX|AVXFMA4|AVX2|AVX512|AVX512MIC|IMCI|KNL|KNC],\
	[Select instructions to be SSE4.0, AVX 1.0, AVX 2.0+FMA, AVX 512, IMCI])],\
	[ac_SIMD=${enable_simd}],[ac_SIMD=GEN])

case ${ax_cv_cxx_compiler_vendor} in
  clang|gnu)
    case ${ac_SIMD} in
      SSE4)
        AC_DEFINE([SSE4],[1],[SSE4 intrinsics])
        SIMD_FLAGS='-msse4.2';;
      AVX)
        AC_DEFINE([AVX1],[1],[AVX intrinsics])
        SIMD_FLAGS='-mavx';;
      AVXFMA4)
        AC_DEFINE([AVXFMA4],[1],[AVX intrinsics with FMA4])
        SIMD_FLAGS='-mavx -mfma4';;
      AVX2)
        AC_DEFINE([AVX2],[1],[AVX2 intrinsics])
        SIMD_FLAGS='-mavx2 -mfma';;
      AVX512|AVX512MIC|KNL)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics])
        SIMD_FLAGS='-mavx512f -mavx512pf -mavx512er -mavx512cd';;
      IMCI|KNC)
        AC_DEFINE([IMCI],[1],[IMCI Intrinsics for Knights Corner])
        SIMD_FLAGS='';;
      GEN)
        AC_DEFINE([GENERIC_VEC],[1],[generic vector code])
        SIMD_FLAGS='';;
      *)
        AC_MSG_ERROR(["SIMD option ${ac_SIMD} not supported by the GCC/Clang compiler"]);;
    esac;;
  intel)
    case ${ac_SIMD} in
      SSE4)
        AC_DEFINE([SSE4],[1],[SSE4 intrinsics])
        SIMD_FLAGS='-msse4.2 -xsse4.2';;
      AVX)
        AC_DEFINE([AVX1],[1],[AVX intrinsics])
        SIMD_FLAGS='-mavx -xavx';;
      AVXFMA4)
        AC_DEFINE([AVXFMA4],[1],[AVX intrinsics with FMA4])
        SIMD_FLAGS='-mavx -xavx -mfma';;
      AVX2)
        AC_DEFINE([AVX2],[1],[AVX2 intrinsics])
        SIMD_FLAGS='-march=core-avx2 -xcore-avx2';;
      AVX512)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics])
        SIMD_FLAGS='-xcore-avx512';;
      AVX512MIC|KNL)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics for Knights Landing])
        SIMD_FLAGS='-xmic-avx512';;
      IMCI|KNC)
        AC_DEFINE([IMCI],[1],[IMCI Intrinsics for Knights Corner])
        SIMD_FLAGS='';;
      GEN)
        AC_DEFINE([GENERIC_VEC],[1],[generic vector code])
        SIMD_FLAGS='';;
      *)
        AC_MSG_ERROR(["SIMD option ${ac_SIMD} not supported by the Intel compiler"]);;
    esac;;
  *)
    AC_MSG_WARN([Compiler unknown, using generic vector code])
    AC_DEFINE([GENERIC_VEC],[1],[generic vector code]);;
esac
AM_CXXFLAGS="$SIMD_FLAGS $AM_CXXFLAGS"
AM_CFLAGS="$SIMD_FLAGS $AM_CFLAGS"

case ${ac_SIMD} in
  AVX512|AVX512MIC|KNL)
    AC_DEFINE([TEST_ZMM],[1],[compile ZMM test]);;
  *)
	;;
esac

############### precision selection
AC_ARG_ENABLE([precision],[AC_HELP_STRING([--enable-precision=single|double],[Select default word size of Real])],[ac_PRECISION=${enable_precision}],[ac_PRECISION=double])
case ${ac_PRECISION} in
     single)
       AC_DEFINE([GRID_DEFAULT_PRECISION_SINGLE],[1],[GRID_DEFAULT_PRECISION is SINGLE] )
     ;;
     double)
       AC_DEFINE([GRID_DEFAULT_PRECISION_DOUBLE],[1],[GRID_DEFAULT_PRECISION is DOUBLE] )
     ;;
esac

############### communication type selection
AC_ARG_ENABLE([comms],[AC_HELP_STRING([--enable-comms=none|mpi|mpi-auto|shmem],[Select communications])],[ac_COMMS=${enable_comms}],[ac_COMMS=none])

case ${ac_COMMS} in
     none)
       AC_DEFINE([GRID_COMMS_NONE],[1],[GRID_COMMS_NONE] )
     ;;
     mpi-auto)
       AC_DEFINE([GRID_COMMS_MPI],[1],[GRID_COMMS_MPI] )
       LX_FIND_MPI
       if test "x$have_CXX_mpi" = 'xno'; then AC_MSG_ERROR(["MPI not found"]); fi
       AM_CXXFLAGS="$MPI_CXXFLAGS $AM_CXXFLAGS"
       AM_CFLAGS="$MPI_CFLAGS $AM_CFLAGS"
       AM_LDFLAGS="`echo $MPI_CXXLDFLAGS | sed -E 's/-l@<:@^ @:>@+//g'` $AM_LDFLAGS"
       LIBS="`echo $MPI_CXXLDFLAGS | sed -E 's/-L@<:@^ @:>@+//g'` $LIBS"
     ;;
     mpi)
       AC_DEFINE([GRID_COMMS_MPI],[1],[GRID_COMMS_MPI] )
     ;;
     shmem)
       AC_DEFINE([GRID_COMMS_SHMEM],[1],[GRID_COMMS_SHMEM] )
     ;;
     *)
     AC_MSG_ERROR([${ac_COMMS} unsupported --enable-comms option]); 
     ;;
esac
AM_CONDITIONAL(BUILD_COMMS_SHMEM,[ test "X${ac_COMMS}X" == "XshmemX" ])
AM_CONDITIONAL(BUILD_COMMS_MPI,[ test "X${ac_COMMS}X" == "XmpiX" || test "X${ac_COMMS}X" == "Xmpi-autoX" ])
AM_CONDITIONAL(BUILD_COMMS_NONE,[ test "X${ac_COMMS}X" == "XnoneX" ])

############### RNG selection
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

############### timer option
AC_ARG_ENABLE([timers],[AC_HELP_STRING([--enable-timers],\
	[Enable system dependent high res timers])],\
	[ac_TIMERS=${enable_timers}],[ac_TIMERS=yes])
case ${ac_TIMERS} in
     yes)
      AC_DEFINE([TIMERS_ON],[1],[TIMERS_ON] )
     ;;
     no)
      AC_DEFINE([TIMERS_OFF],[1],[TIMERS_OFF] )
     ;;
     *)
      AC_MSG_ERROR([${ac_TIMERS} unsupported --enable-timers option]); 
     ;;
esac

############### Chroma regression test
AC_ARG_ENABLE([chroma],[AC_HELP_STRING([--enable-chroma],[Expect chroma compiled under c++11 ])],ac_CHROMA=yes,ac_CHROMA=no)
case ${ac_CHROMA} in
     yes|no)
     ;;
     *)
       AC_MSG_ERROR([${ac_CHROMA} unsupported --enable-chroma option]); 
     ;;
esac
AM_CONDITIONAL(BUILD_CHROMA_REGRESSION,[ test "X${ac_CHROMA}X" == "XyesX" ])

############### Doxygen
AC_PROG_DOXYGEN

if test -n "$DOXYGEN"
then
AC_CONFIG_FILES([docs/doxy.cfg])
fi

############### Ouput
cwd=`pwd -P`; cd ${srcdir}; abs_srcdir=`pwd -P`; cd ${cwd}
AM_CXXFLAGS="-I${abs_srcdir}/include $AM_CXXFLAGS"
AM_CFLAGS="-I${abs_srcdir}/include $AM_CFLAGS"
AM_LDFLAGS="-L${cwd}/lib $AM_LDFLAGS"
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(tests/Makefile)
AC_CONFIG_FILES(tests/IO/Makefile)
AC_CONFIG_FILES(tests/core/Makefile)
AC_CONFIG_FILES(tests/debug/Makefile)
AC_CONFIG_FILES(tests/forces/Makefile)
AC_CONFIG_FILES(tests/hmc/Makefile)
AC_CONFIG_FILES(tests/solver/Makefile)
AC_CONFIG_FILES(tests/qdpxx/Makefile)
AC_CONFIG_FILES(benchmarks/Makefile)
AC_OUTPUT

echo "
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary of configuration for $PACKAGE v$VERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

----- PLATFORM ----------------------------------------
- architecture (build)          : $build_cpu
- os (build)                    : $build_os
- architecture (target)         : $target_cpu
- os (target)                   : $target_os
- compiler vendor               : ${ax_cv_cxx_compiler_vendor}
- compiler version              : ${ax_cv_gxx_version}
----- BUILD OPTIONS -----------------------------------
- SIMD                          : ${ac_SIMD}
- communications type           : ${ac_COMMS}
- default precision             : ${ac_PRECISION}
- RNG choice                    : ${ac_RNG} 
- GMP                           : `if test "x$have_gmp" = xtrue; then echo yes; else echo no; fi`
- LAPACK                        : ${ac_LAPACK}
- build DOXYGEN documentation   : `if test "x$enable_doc" = xyes; then echo yes; else echo no; fi`
- graphs and diagrams           : `if test "x$enable_dot" = xyes; then echo yes; else echo no; fi`
----- BUILD FLAGS -------------------------------------
- CXXFLAGS: 
`echo ${AM_CXXFLAGS} ${CXXFLAGS} | sed 's/ -/\n\t-/g' | sed 's/^-/\t-/g'`
- LDFLAGS:
`echo ${AM_LDFLAGS} ${LDFLAGS} | sed 's/ -/\n\t-/g' | sed 's/^-/\t-/g'`
- LIBS:
`echo ${LIBS} | sed 's/ -/\n\t-/g' | sed 's/^-/\t-/g'`
-------------------------------------------------------
"
