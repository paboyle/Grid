AC_PREREQ([2.63])
AC_INIT([Grid], [0.7.0], [https://github.com/paboyle/Grid], [Grid])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([subdir-objects 1.13])
AM_EXTRA_RECURSIVE_TARGETS([tests bench])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([lib/Grid.h])
AC_CONFIG_HEADERS([lib/Config.h],[sed -i 's|PACKAGE_|GRID_|' lib/Config.h])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

################ Get git info
#AC_REVISION([m4_esyscmd_s([./scripts/configure.commit])])

############### Checks for programs
AC_PROG_CXX
AC_PROG_RANLIB

############### Get compiler informations
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
AX_COMPILER_VENDOR
AC_DEFINE_UNQUOTED([CXX_COMP_VENDOR],["$ax_cv_cxx_compiler_vendor"],
      [vendor of C++ compiler that will compile the code])
AX_GXX_VERSION
AC_DEFINE_UNQUOTED([GXX_VERSION],["$GXX_VERSION"],
      [version of g++ that will compile the code])

CXXFLAGS="-O3 $CXXFLAGS"


############### Checks for typedefs, structures, and compiler characteristics
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

############### OpenMP
AC_OPENMP
ac_openmp=no
if test "${OPENMP_CXXFLAGS}X" != "X"; then
  ac_openmp=yes
  AM_CXXFLAGS="$OPENMP_CXXFLAGS $AM_CXXFLAGS"
  AM_LDFLAGS="$OPENMP_CXXFLAGS $AM_LDFLAGS"
fi

############### Checks for header files
AC_CHECK_HEADERS(stdint.h)
AC_CHECK_HEADERS(mm_malloc.h)
AC_CHECK_HEADERS(malloc/malloc.h)
AC_CHECK_HEADERS(malloc.h)
AC_CHECK_HEADERS(endian.h)
AC_CHECK_HEADERS(execinfo.h)
AC_CHECK_DECLS([ntohll],[], [], [[#include <arpa/inet.h>]])
AC_CHECK_DECLS([be64toh],[], [], [[#include <arpa/inet.h>]])

############### GMP and MPFR
AC_ARG_WITH([gmp],
    [AS_HELP_STRING([--with-gmp=prefix],
    [try this for a non-standard install prefix of the GMP library])],
    [AM_CXXFLAGS="-I$with_gmp/include $AM_CXXFLAGS"]
    [AM_LDFLAGS="-L$with_gmp/lib $AM_LDFLAGS"])
AC_ARG_WITH([mpfr],
    [AS_HELP_STRING([--with-mpfr=prefix],
    [try this for a non-standard install prefix of the MPFR library])],
    [AM_CXXFLAGS="-I$with_mpfr/include $AM_CXXFLAGS"]
    [AM_LDFLAGS="-L$with_mpfr/lib $AM_LDFLAGS"])

############### FFTW3
AC_ARG_WITH([fftw],
            [AS_HELP_STRING([--with-fftw=prefix],
            [try this for a non-standard install prefix of the FFTW3 library])],
            [AM_CXXFLAGS="-I$with_fftw/include $AM_CXXFLAGS"]
            [AM_LDFLAGS="-L$with_fftw/lib $AM_LDFLAGS"])

############### LIME
AC_ARG_WITH([lime],
            [AS_HELP_STRING([--with-lime=prefix],
            [try this for a non-standard install prefix of the LIME library])],
            [AM_CXXFLAGS="-I$with_lime/include $AM_CXXFLAGS"]
            [AM_LDFLAGS="-L$with_lime/lib $AM_LDFLAGS"])

############### lapack
AC_ARG_ENABLE([lapack],
    [AC_HELP_STRING([--enable-lapack=yes|no|prefix], [enable LAPACK])],
    [ac_LAPACK=${enable_lapack}], [ac_LAPACK=no])

case ${ac_LAPACK} in
    no)
        ;;
    yes)
        AC_DEFINE([USE_LAPACK],[1],[use LAPACK]);;
    *)
        AM_CXXFLAGS="-I$ac_LAPACK/include $AM_CXXFLAGS"
        AM_LDFLAGS="-L$ac_LAPACK/lib $AM_LDFLAGS"
        AC_DEFINE([USE_LAPACK],[1],[use LAPACK]);;
esac

############### FP16 conversions
AC_ARG_ENABLE([sfw-fp16],
    [AC_HELP_STRING([--enable-sfw-fp16=yes|no], [enable software fp16 comms])],
    [ac_SFW_FP16=${enable_sfw_fp16}], [ac_SFW_FP16=yes])
case ${ac_SFW_FP16} in
    yes)
      AC_DEFINE([SFW_FP16],[1],[software conversion to fp16]);;
    no);;
    *)
      AC_MSG_ERROR(["SFW FP16 option not supported ${ac_SFW_FP16}"]);;
esac

############### MKL
AC_ARG_ENABLE([mkl],
    [AC_HELP_STRING([--enable-mkl=yes|no|prefix], [enable Intel MKL for LAPACK & FFTW])],
    [ac_MKL=${enable_mkl}], [ac_MKL=no])

case ${ac_MKL} in
    no)
        ;;
    yes)
        AC_DEFINE([USE_MKL], [1], [Define to 1 if you use the Intel MKL]);;
    *)
        AM_CXXFLAGS="-I$ac_MKL/include $AM_CXXFLAGS"
        AM_LDFLAGS="-L$ac_MKL/lib $AM_LDFLAGS"
        AC_DEFINE([USE_MKL], [1], [Define to 1 if you use the Intel MKL]);;
esac

############### HDF5
AC_ARG_WITH([hdf5],
    [AS_HELP_STRING([--with-hdf5=prefix],
    [try this for a non-standard install prefix of the HDF5 library])],
    [AM_CXXFLAGS="-I$with_hdf5/include $AM_CXXFLAGS"]
    [AM_LDFLAGS="-L$with_hdf5/lib $AM_LDFLAGS"])

############### first-touch
AC_ARG_ENABLE([numa],
    [AC_HELP_STRING([--enable-numa=yes|no|prefix], [enable first touch numa opt])],
    [ac_NUMA=${enable_NUMA}],[ac_NUMA=no])

case ${ac_NUMA} in
    no)
        ;;
    yes)
        AC_DEFINE([GRID_NUMA],[1],[First touch numa locality]);;
    *)
        AC_DEFINE([GRID_NUMA],[1],[First touch numa locality]);;
esac

############### Checks for library functions
CXXFLAGS_CPY=$CXXFLAGS
LDFLAGS_CPY=$LDFLAGS
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"

AC_CHECK_FUNCS([gettimeofday])

if test "${ac_MKL}x" != "nox"; then
    AC_SEARCH_LIBS([mkl_set_interface_layer], [mkl_rt], [],
                   [AC_MSG_ERROR("MKL enabled but library not found")])
fi

AC_SEARCH_LIBS([__gmpf_init], [gmp],
               [AC_SEARCH_LIBS([mpfr_init], [mpfr],
                               [AC_DEFINE([HAVE_LIBMPFR], [1],
                                          [Define to 1 if you have the `MPFR' library])]
                               [have_mpfr=true], [AC_MSG_ERROR([MPFR library not found])])]
               [AC_DEFINE([HAVE_LIBGMP], [1], [Define to 1 if you have the `GMP' library])]
               [have_gmp=true])

if test "${ac_LAPACK}x" != "nox"; then
    AC_SEARCH_LIBS([LAPACKE_sbdsdc], [lapack], [],
                   [AC_MSG_ERROR("LAPACK enabled but library not found")])
fi

AC_SEARCH_LIBS([fftw_execute], [fftw3],
               [AC_SEARCH_LIBS([fftwf_execute], [fftw3f], [],
                               [AC_MSG_ERROR("single precision FFTW library not found")])]
               [AC_DEFINE([HAVE_FFTW], [1], [Define to 1 if you have the `FFTW' library])]
               [have_fftw=true])

AC_SEARCH_LIBS([limeCreateReader], [lime],
               [AC_DEFINE([HAVE_LIME], [1], [Define to 1 if you have the `LIME' library])]
               [have_lime=true],
	       [AC_MSG_WARN(C-LIME library was not found in your system.
In order to use ILGG file format please install or provide the correct path to your installation
Info at: http://usqcd.jlab.org/usqcd-docs/c-lime/)])


AC_SEARCH_LIBS([H5Fopen], [hdf5_cpp],
               [AC_DEFINE([HAVE_HDF5], [1], [Define to 1 if you have the `HDF5' library])]
               [have_hdf5=true]
               [LIBS="${LIBS} -lhdf5"], [], [-lhdf5])
AM_CONDITIONAL(BUILD_HDF5, [ test "${have_hdf5}X" == "trueX" ])

CXXFLAGS=$CXXFLAGS_CPY
LDFLAGS=$LDFLAGS_CPY

############### SIMD instruction selection
AC_ARG_ENABLE([simd],[AC_HELP_STRING([--enable-simd=code],
              [select SIMD target (cf. README.md)])], [ac_SIMD=${enable_simd}], [ac_SIMD=GEN])

AC_ARG_ENABLE([gen-simd-width],
            [AS_HELP_STRING([--enable-gen-simd-width=size],
            [size (in bytes) of the generic SIMD vectors (default: 32)])],
            [ac_gen_simd_width=$enable_gen_simd_width],
            [ac_gen_simd_width=32])

case ${ax_cv_cxx_compiler_vendor} in
  clang|gnu)
    case ${ac_SIMD} in
      SSE4)
        AC_DEFINE([SSE4],[1],[SSE4 intrinsics])
	case ${ac_SFW_FP16} in
	  yes)
	  SIMD_FLAGS='-msse4.2';;
	  no)
	  SIMD_FLAGS='-msse4.2 -mf16c';;
	  *)
          AC_MSG_ERROR(["SFW_FP16 must be either yes or no value ${ac_SFW_FP16} "]);;
	esac;;
      AVX)
        AC_DEFINE([AVX1],[1],[AVX intrinsics])
        SIMD_FLAGS='-mavx -mf16c';;
      AVXFMA4)
        AC_DEFINE([AVXFMA4],[1],[AVX intrinsics with FMA4])
        SIMD_FLAGS='-mavx -mfma4 -mf16c';;
      AVXFMA)
        AC_DEFINE([AVXFMA],[1],[AVX intrinsics with FMA3])
        SIMD_FLAGS='-mavx -mfma -mf16c';;
      AVX2)
        AC_DEFINE([AVX2],[1],[AVX2 intrinsics])
        SIMD_FLAGS='-mavx2 -mfma -mf16c';;
      AVX512)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics])
        SIMD_FLAGS='-mavx512f -mavx512pf -mavx512er -mavx512cd';;
      KNC)
        AC_DEFINE([IMCI],[1],[IMCI intrinsics for Knights Corner])
        SIMD_FLAGS='';;
      KNL)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics])
        SIMD_FLAGS='-march=knl';;
      GEN)
        AC_DEFINE([GEN],[1],[generic vector code])
        AC_DEFINE_UNQUOTED([GEN_SIMD_WIDTH],[$ac_gen_simd_width],
                           [generic SIMD vector width (in bytes)])
        SIMD_GEN_WIDTH_MSG=" (width= $ac_gen_simd_width)"
        SIMD_FLAGS='';;
      NEONv8)
        AC_DEFINE([NEONV8],[1],[ARMv8 NEON])
        SIMD_FLAGS='';;
      QPX|BGQ)
        AC_DEFINE([QPX],[1],[QPX intrinsics for BG/Q])
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
      AVXFMA)
        AC_DEFINE([AVXFMA],[1],[AVX intrinsics with FMA3])
        SIMD_FLAGS='-mavx -fma';;
      AVX2)
        AC_DEFINE([AVX2],[1],[AVX2 intrinsics])
        SIMD_FLAGS='-march=core-avx2 -xcore-avx2';;
      AVX512)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics])
        SIMD_FLAGS='-xcore-avx512';;
      KNC)
        AC_DEFINE([IMCI],[1],[IMCI Intrinsics for Knights Corner])
        SIMD_FLAGS='';;
      KNL)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics for Knights Landing])
        SIMD_FLAGS='-xmic-avx512';;
      GEN)
        AC_DEFINE([GEN],[1],[generic vector code])
        AC_DEFINE_UNQUOTED([GEN_SIMD_WIDTH],[$ac_gen_simd_width],
                           [generic SIMD vector width (in bytes)])
        SIMD_GEN_WIDTH_MSG=" (width= $ac_gen_simd_width)"
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
  AVX512|KNL)
    AC_DEFINE([TEST_ZMM],[1],[compile ZMM test]);;
  *)
	;;
esac

############### Precision selection
AC_ARG_ENABLE([precision],
              [AC_HELP_STRING([--enable-precision=single|double],
                              [Select default word size of Real])],
              [ac_PRECISION=${enable_precision}],[ac_PRECISION=double])

case ${ac_PRECISION} in
     single)
       AC_DEFINE([GRID_DEFAULT_PRECISION_SINGLE],[1],[GRID_DEFAULT_PRECISION is SINGLE] )
     ;;
     double)
       AC_DEFINE([GRID_DEFAULT_PRECISION_DOUBLE],[1],[GRID_DEFAULT_PRECISION is DOUBLE] )
     ;;
esac

############### communication type selection
AC_ARG_ENABLE([comms],[AC_HELP_STRING([--enable-comms=none|mpi|mpi-auto|mpi3|mpi3-auto|shmem],
              [Select communications])],[ac_COMMS=${enable_comms}],[ac_COMMS=none])

case ${ac_COMMS} in
     none)
        AC_DEFINE([GRID_COMMS_NONE],[1],[GRID_COMMS_NONE] )
        comms_type='none'
     ;;
     mpi3l*)
       AC_DEFINE([GRID_COMMS_MPI3L],[1],[GRID_COMMS_MPI3L] )
       comms_type='mpi3l'
     ;;
     mpi3*)
        AC_DEFINE([GRID_COMMS_MPI3],[1],[GRID_COMMS_MPI3] )
        comms_type='mpi3'
     ;;
     mpi*)
        AC_DEFINE([GRID_COMMS_MPI],[1],[GRID_COMMS_MPI] )
        comms_type='mpi'
     ;;
     shmem)
        AC_DEFINE([GRID_COMMS_SHMEM],[1],[GRID_COMMS_SHMEM] )
        comms_type='shmem'
     ;;
     *)
        AC_MSG_ERROR([${ac_COMMS} unsupported --enable-comms option]);
     ;;
esac
case ${ac_COMMS} in
    *-auto)
        LX_FIND_MPI
        if test "x$have_CXX_mpi" = 'xno'; then AC_MSG_ERROR(["The configure could not find the MPI compilation flags. N.B. The -auto mode is not supported by Cray wrappers. Use the non -auto version in this case."]); fi
        AM_CXXFLAGS="$MPI_CXXFLAGS $AM_CXXFLAGS"
        AM_CFLAGS="$MPI_CFLAGS $AM_CFLAGS"
        AM_LDFLAGS="`echo $MPI_CXXLDFLAGS | sed -E 's/-l@<:@^ @:>@+//g'` $AM_LDFLAGS"
        LIBS="`echo $MPI_CXXLDFLAGS | sed -E 's/-L@<:@^ @:>@+//g'` $LIBS";;
    *)
        ;;
esac

AM_CONDITIONAL(BUILD_COMMS_SHMEM, [ test "${comms_type}X" == "shmemX" ])
AM_CONDITIONAL(BUILD_COMMS_MPI,   [ test "${comms_type}X" == "mpiX" ])
AM_CONDITIONAL(BUILD_COMMS_MPI3,  [ test "${comms_type}X" == "mpi3X" ] )
AM_CONDITIONAL(BUILD_COMMS_MPI3L, [ test "${comms_type}X" == "mpi3lX" ] )
AM_CONDITIONAL(BUILD_COMMS_NONE,  [ test "${comms_type}X" == "noneX" ])

############### RNG selection
AC_ARG_ENABLE([rng],[AC_HELP_STRING([--enable-rng=ranlux48|mt19937|sitmo],\
	            [Select Random Number Generator to be used])],\
	            [ac_RNG=${enable_rng}],[ac_RNG=sitmo])

case ${ac_RNG} in
     ranlux48)
      AC_DEFINE([RNG_RANLUX],[1],[RNG_RANLUX] )
     ;;
     mt19937)
      AC_DEFINE([RNG_MT19937],[1],[RNG_MT19937] )
     ;;
     sitmo)
      AC_DEFINE([RNG_SITMO],[1],[RNG_SITMO] )
     ;;
     *)
      AC_MSG_ERROR([${ac_RNG} unsupported --enable-rng option]);
     ;;
esac

############### Timer option
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
AC_ARG_ENABLE([chroma],[AC_HELP_STRING([--enable-chroma],
              [Expect chroma compiled under c++11 ])],ac_CHROMA=yes,ac_CHROMA=no)

case ${ac_CHROMA} in
     yes|no)
     ;;
     *)
       AC_MSG_ERROR([${ac_CHROMA} unsupported --enable-chroma option]);
     ;;
esac

AM_CONDITIONAL(BUILD_CHROMA_REGRESSION,[ test "X${ac_CHROMA}X" == "XyesX" ])

############### Doxygen
DX_DOXYGEN_FEATURE([OFF])
DX_DOT_FEATURE([OFF])
DX_HTML_FEATURE([ON])
DX_CHM_FEATURE([OFF])
DX_CHI_FEATURE([OFF])
DX_MAN_FEATURE([OFF])
DX_RTF_FEATURE([OFF])
DX_XML_FEATURE([OFF])
DX_PDF_FEATURE([OFF])
DX_PS_FEATURE([OFF])
DX_INIT_DOXYGEN([$PACKAGE_NAME], [doxygen.cfg])

############### Ouput
cwd=`pwd -P`; cd ${srcdir}; abs_srcdir=`pwd -P`; cd ${cwd}
GRID_CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
GRID_LDFLAGS="$AM_LDFLAGS $LDFLAGS"
GRID_LIBS=$LIBS
GRID_SHORT_SHA=`git rev-parse --short HEAD`
GRID_SHA=`git rev-parse HEAD`
GRID_BRANCH=`git rev-parse --abbrev-ref HEAD`
AM_CXXFLAGS="-I${abs_srcdir}/include $AM_CXXFLAGS"
AM_CFLAGS="-I${abs_srcdir}/include $AM_CFLAGS"
AM_LDFLAGS="-L${cwd}/lib $AM_LDFLAGS"
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([GRID_CXXFLAGS])
AC_SUBST([GRID_LDFLAGS])
AC_SUBST([GRID_LIBS])
AC_SUBST([GRID_SHA])
AC_SUBST([GRID_BRANCH])

git_commit=`cd $srcdir && ./scripts/configure.commit`

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary of configuration for $PACKAGE v$VERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
----- GIT VERSION -------------------------------------
$git_commit
----- PLATFORM ----------------------------------------
architecture (build)        : $build_cpu
os (build)                  : $build_os
architecture (target)       : $target_cpu
os (target)                 : $target_os
compiler vendor             : ${ax_cv_cxx_compiler_vendor}
compiler version            : ${ax_cv_gxx_version}
----- BUILD OPTIONS -----------------------------------
SIMD                        : ${ac_SIMD}${SIMD_GEN_WIDTH_MSG}
Threading                   : ${ac_openmp}
Communications type         : ${comms_type}
Default precision           : ${ac_PRECISION}
Software FP16 conversion    : ${ac_SFW_FP16}
RNG choice                  : ${ac_RNG}
GMP                         : `if test "x$have_gmp" = xtrue; then echo yes; else echo no; fi`
LAPACK                      : ${ac_LAPACK}
FFTW                        : `if test "x$have_fftw" = xtrue; then echo yes; else echo no; fi`
LIME (ILDG support)         : `if test "x$have_lime" = xtrue; then echo yes; else echo no; fi`
HDF5                        : `if test "x$have_hdf5" = xtrue; then echo yes; else echo no; fi`
build DOXYGEN documentation : `if test "$DX_FLAG_doc" = '1'; then echo yes; else echo no; fi`
----- BUILD FLAGS -------------------------------------
CXXFLAGS:
`echo ${AM_CXXFLAGS} ${CXXFLAGS} | tr ' ' '\n' | sed 's/^-/    -/g'`
LDFLAGS:
`echo ${AM_LDFLAGS} ${LDFLAGS} | tr ' ' '\n' | sed 's/^-/    -/g'`
LIBS:
`echo ${LIBS} | tr ' ' '\n' | sed 's/^-/    -/g'`
-------------------------------------------------------" > grid.configure.summary

GRID_SUMMARY="`cat grid.configure.summary`"
AM_SUBST_NOTMAKE([GRID_SUMMARY])
AC_SUBST([GRID_SUMMARY])

AC_CONFIG_FILES([grid-config], [chmod +x grid-config])
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(tests/Makefile)
AC_CONFIG_FILES(tests/IO/Makefile)
AC_CONFIG_FILES(tests/core/Makefile)
AC_CONFIG_FILES(tests/debug/Makefile)
AC_CONFIG_FILES(tests/forces/Makefile)
AC_CONFIG_FILES(tests/hadrons/Makefile)
AC_CONFIG_FILES(tests/hmc/Makefile)
AC_CONFIG_FILES(tests/solver/Makefile)
AC_CONFIG_FILES(tests/smearing/Makefile)
AC_CONFIG_FILES(tests/qdpxx/Makefile)
AC_CONFIG_FILES(tests/testu01/Makefile)
AC_CONFIG_FILES(benchmarks/Makefile)
AC_CONFIG_FILES(extras/Makefile)
AC_CONFIG_FILES(extras/Hadrons/Makefile)
AC_OUTPUT

echo ""
cat grid.configure.summary
echo ""

