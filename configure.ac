AC_PREREQ([2.69])
AC_INIT([Grid],[0.7.0],[https://github.com/paboyle/Grid],[Grid])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([subdir-objects 1.13])
AM_EXTRA_RECURSIVE_TARGETS([tests bench])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([Grid/Grid.h])
AC_CONFIG_HEADERS([Grid/Config.h],[[$SED_INPLACE -e 's|PACKAGE_|GRID_|' -e 's|[[:space:]]PACKAGE[[:space:]]| GRID_PACKAGE |' -e 's|[[:space:]]VERSION[[:space:]]| GRID_PACKAGE_VERSION |' Grid/Config.h]],
    [if test x"$host_os" == x"${host_os#darwin}" ; then]
        [SED_INPLACE="sed -i"]
    [else]
        [SED_INPLACE="sed -i .bak"]
    [fi])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

################ Get git info
#AC_REVISION([m4_esyscmd_s([./scripts/configure.commit])])

################ Set flags
# do not move!
AC_ARG_ENABLE([debug],[AS_HELP_STRING([--enable-debug=yes|no],[enable debug compilation ])],
    [ac_DEBUG=${enable_debug}], [ac_DEBUG=no])
case ${ac_DEBUG} in
    yes)
    CXXFLAGS="-g $CXXFLAGS "
        ;;
    no)
    CXXFLAGS="-O3 $CXXFLAGS "
        ;;
    *)
    CXXFLAGS="-O3 $CXXFLAGS "
        ;;
esac


############### Checks for programs
AC_PROG_CXX
AC_PROG_RANLIB

############### Get compiler informations
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX(17,noext,mandatory)
AX_COMPILER_VENDOR
AC_DEFINE_UNQUOTED([CXX_COMP_VENDOR],["$ax_cv_cxx_compiler_vendor"],
      [vendor of C++ compiler that will compile the code])
AX_GXX_VERSION
AC_DEFINE_UNQUOTED([GXX_VERSION],["$GXX_VERSION"],
      [version of g++ that will compile the code])


############### Checks for typedefs, structures, and compiler characteristics
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

############### OpenMP
AC_OPENMP
ac_openmp=no
if test "${ac_cv_prog_cxx_openmp}X" != "noX"; then
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

############## Standard libraries
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([stdc++],[abort])

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

############### OpenSSL
AC_ARG_WITH([openssl],
            [AS_HELP_STRING([--with-openssl=prefix],
            [try this for a non-standard install prefix of the OpenSSL library])],
            [AM_CXXFLAGS="-I$with_openssl/include $AM_CXXFLAGS"]
            [AM_LDFLAGS="-L$with_openssl/lib $AM_LDFLAGS"])

############### lapack
AC_ARG_ENABLE([lapack],
    [AS_HELP_STRING([--enable-lapack=yes|no|prefix],[enable LAPACK])],
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

############### tracing
AC_ARG_ENABLE([tracing],
    [AS_HELP_STRING([--enable-tracing=none|nvtx|roctx|timer],[enable tracing])],
    [ac_TRACING=${enable_tracing}], [ac_TRACING=none])

case ${ac_TRACING} in
    nvtx)
        AC_DEFINE([GRID_TRACING_NVTX],[1],[use NVTX])
	LIBS="${LIBS} -lnvToolsExt"
	;;
    roctx)
        AC_DEFINE([GRID_TRACING_ROCTX],[1],[use ROCTX])
	LIBS="${LIBS} -lroctx64"
	;;
    timer)
        AC_DEFINE([GRID_TRACING_TIMER],[1],[use TIMER]);;
    *)
	AC_DEFINE([GRID_TRACING_NONE],[1],[no tracing]);;
esac

############### fermions
AC_ARG_ENABLE([fermion-reps],
     [AS_HELP_STRING([--enable-fermion-reps=yes|no],[enable extra fermion representation support])],
     [ac_FERMION_REPS=${enable_fermion_reps}], [ac_FERMION_REPS=yes])

AM_CONDITIONAL(BUILD_FERMION_REPS, [ test "${ac_FERMION_REPS}X" == "yesX" ])

AC_ARG_ENABLE([gparity],
     [AS_HELP_STRING([--enable-gparity=yes|no],[enable G-parity support])],
     [ac_GPARITY=${enable_gparity}], [ac_GPARITY=yes])

AM_CONDITIONAL(BUILD_GPARITY, [ test "${ac_GPARITY}X" == "yesX" ])

AC_ARG_ENABLE([zmobius],
     [AS_HELP_STRING([--enable-zmobius=yes|no],[enable Zmobius support])],
     [ac_ZMOBIUS=${enable_zmobius}], [ac_ZMOBIUS=yes])

AM_CONDITIONAL(BUILD_ZMOBIUS, [ test "${ac_ZMOBIUS}X" == "yesX" ])


case ${ac_FERMION_REPS} in
   yes) AC_DEFINE([ENABLE_FERMION_REPS],[1],[non QCD fermion reps]);;
esac
case ${ac_GPARITY} in
   yes) AC_DEFINE([ENABLE_GPARITY],[1],[fermion actions with GPARITY BCs]);;
esac
case ${ac_ZMOBIUS} in
   yes) AC_DEFINE([ENABLE_ZMOBIUS],[1],[Zmobius fermion actions]);;
esac
############### Nc
AC_ARG_ENABLE([Nc],
    [AS_HELP_STRING([--enable-Nc=2|3|4|5],[enable number of colours])],
    [ac_Nc=${enable_Nc}], [ac_Nc=3])

case ${ac_Nc} in
    2)
        AC_DEFINE([Config_Nc],[2],[Gauge group Nc]);;
    3)
        AC_DEFINE([Config_Nc],[3],[Gauge group Nc]);;
    4)
        AC_DEFINE([Config_Nc],[4],[Gauge group Nc]);;
    5)
        AC_DEFINE([Config_Nc],[5],[Gauge group Nc]);;
    8)
        AC_DEFINE([Config_Nc],[8],[Gauge group Nc]);;
    *)
      AC_MSG_ERROR(["Unsupport gauge group choice Nc = ${ac_Nc}"]);;
esac

############### Symplectic group
AC_ARG_ENABLE([Sp],
    [AC_HELP_STRING([--enable-Sp=yes|no], [enable gauge group Sp2n])],
    [ac_ENABLE_SP=${enable_Sp}], [ac_ENABLE_SP=no])

AM_CONDITIONAL(BUILD_SP, [ test "${ac_ENABLE_SP}X" == "yesX" ])

case ${ac_ENABLE_SP} in
   yes)
        AC_DEFINE([Sp2n_config],[1],[gauge group Sp2n], [have_sp2n=true]);;
   no)
        AC_DEFINE([Sp2n_config],[0],[gauge group SUn], [have_sp2n=false]);;
    *)
        AC_MSG_ERROR(["--enable-Sp is either yes or no"]);;
esac

############### FP16 conversions
AC_ARG_ENABLE([sfw-fp16],
    [AS_HELP_STRING([--enable-sfw-fp16=yes|no],[enable software fp16 comms])],
    [ac_SFW_FP16=${enable_sfw_fp16}], [ac_SFW_FP16=yes])
case ${ac_SFW_FP16} in
    yes)
      AC_DEFINE([SFW_FP16],[1],[software conversion to fp16]);;
    no);;
    *)
      AC_MSG_ERROR(["SFW FP16 option not supported ${ac_SFW_FP16}"]);;
esac

############### Default to accelerator cshift, but revert to host if UCX is buggy or other reasons
AC_ARG_ENABLE([accelerator-aware-mpi],
    [AS_HELP_STRING([--enable-accelerator-aware-mpi=yes|no],[run mpi transfers from device])],
    [ac_ACCELERATOR_AWARE_MPI=${enable_accelerator_aware_mpi}], [ac_ACCELERATOR_AWARE_MPI=yes])

case ${ac_ACCELERATOR_AWARE_MPI} in
    yes)
      AC_DEFINE([ACCELERATOR_CSHIFT],[1],[ Cshift runs on host])
      AC_DEFINE([ACCELERATOR_AWARE_MPI],[1],[ Stencil can use device pointers]);;
    *);;
esac


############### SYCL/CUDA/HIP/none
AC_ARG_ENABLE([accelerator],
    [AS_HELP_STRING([--enable-accelerator=cuda|sycl|hip|none],[enable none,cuda,sycl,hip acceleration])],
    [ac_ACCELERATOR=${enable_accelerator}], [ac_ACCELERATOR=none])
case ${ac_ACCELERATOR} in
    cuda)
      echo CUDA acceleration
      LIBS="${LIBS} -lcuda"
      AC_DEFINE([GRID_CUDA],[1],[Use CUDA offload]);;
    sycl)
      echo SYCL acceleration
      AC_DEFINE([GRID_SYCL],[1],[Use SYCL offload]);;
    hip)
      echo HIP acceleration
      AC_DEFINE([GRID_HIP],[1],[Use HIP offload]);;
    none)
      echo NO acceleration    ;;
    no)
      echo NO acceleration    ;;
    *)
      AC_MSG_ERROR(["Acceleration not suppoorted ${ac_ACCELERATOR}"]);;
esac

############### UNIFIED MEMORY
AC_ARG_ENABLE([unified],
    [AS_HELP_STRING([--enable-unified=yes|no],[enable unified address space for accelerator loops])],
    [ac_UNIFIED=${enable_unified}], [ac_UNIFIED=yes])
case ${ac_UNIFIED} in
    yes)
      echo Unified memory for accelerator loops
      AC_DEFINE([GRID_UVM],[1],[Use unified address space]);;
    no)
      echo Manual memory copy for accelerator loops;;
    *)
      AC_MSG_ERROR(["Unified virtual memory option not suppoorted ${ac_UNIFIED}"]);;
esac

############### Intel libraries
AC_ARG_ENABLE([mkl],
    [AS_HELP_STRING([--enable-mkl=yes|no|prefix],[enable Intel MKL for LAPACK & FFTW])],
    [ac_MKL=${enable_mkl}], [ac_MKL=no])
AC_ARG_ENABLE([ipp],
    [AS_HELP_STRING([--enable-ipp=yes|no|prefix],[enable Intel IPP for fast CRC32C])],
    [ac_IPP=${enable_ipp}], [ac_IPP=no])

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

case ${ac_IPP} in
    no)
        ;;
    yes)
        AC_DEFINE([USE_IPP], [1], [Define to 1 if you use the Intel IPP]);;
    *)
        AM_CXXFLAGS="-I$ac_IPP/include $AM_CXXFLAGS"
        AM_LDFLAGS="-L$ac_IPP/lib $AM_LDFLAGS"
        AC_DEFINE([USE_IPP], [1], [Define to 1 if you use the Intel IPP]);;
esac

############### HDF5
AC_ARG_WITH([hdf5],
    [AS_HELP_STRING([--with-hdf5=prefix],
    [try this for a non-standard install prefix of the HDF5 library])],
    [AM_CXXFLAGS="-I$with_hdf5/include $AM_CXXFLAGS"]
    [AM_LDFLAGS="-L$with_hdf5/lib $AM_LDFLAGS"])


############### Checks for library functions
CXXFLAGS_CPY=$CXXFLAGS
LDFLAGS_CPY=$LDFLAGS
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"

AC_CHECK_FUNCS([gettimeofday])

if test "${ac_MKL}x" != "nox"; then
    AC_SEARCH_LIBS([mkl_set_interface_layer], [mkl_rt], [],
                   [AC_MSG_ERROR("Intel MKL enabled but library not found")])
fi

if test "${ac_IPP}x" != "nox"; then
     AC_SEARCH_LIBS([ippsCRC32C_8u], [ippdc],
                   [LIBS="${LIBS} -lippdc -lippvm -lipps -lippcore"],
                   [AC_MSG_ERROR("Intel IPP enabled but library not found")])
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
	             [AC_MSG_WARN(LIME library was not found in your system.)])

AC_SEARCH_LIBS([SHA256_Init], [crypto],
               [AC_DEFINE([HAVE_CRYPTO], [1], [Define to 1 if you have the `OpenSSL' library])]
               [have_crypto=true],
	             [AC_MSG_ERROR(OpenSSL library was not found in your system.)])
AC_CHECK_HEADER([openssl/sha.h], [], [AC_MSG_ERROR(OpenSSL library found but without headers.)], [AC_INCLUDES_DEFAULT([])])

AC_SEARCH_LIBS([crc32], [z],
               [AC_DEFINE([HAVE_ZLIB], [1], [Define to 1 if you have the `LIBZ' library])]
               [have_zlib=true] [LIBS="${LIBS} -lz"],
	       [AC_MSG_ERROR(zlib library was not found in your system.)])


AC_SEARCH_LIBS([H5Fopen], [hdf5_cpp],
               [AC_DEFINE([HAVE_HDF5], [1], [Define to 1 if you have the `HDF5' library])]
               [have_hdf5=true]
               [LIBS="${LIBS} -lhdf5"], [], [-lhdf5])
AM_CONDITIONAL(BUILD_HDF5, [ test "${have_hdf5}X" == "trueX" ])

CXXFLAGS=$CXXFLAGS_CPY
LDFLAGS=$LDFLAGS_CPY

############### SIMD instruction selection
AC_ARG_ENABLE([simd],[AS_HELP_STRING([--enable-simd=code],[select SIMD target (cf. README.md)])], [ac_SIMD=${enable_simd}], [ac_SIMD=GEN])

AC_ARG_ENABLE([gen-simd-width],
            [AS_HELP_STRING([--enable-gen-simd-width=size],
            [size (in bytes) of the generic SIMD vectors (default: 64)])],
            [ac_gen_simd_width=$enable_gen_simd_width],
            [ac_gen_simd_width=64])

AC_ARG_ENABLE([gen-scalar],
            [AS_HELP_STRING([--enable-gen-scalar=yes|no],
            [enable generic scalar implementation])],
            [ac_gen_scalar=$enable_gen_scalar],
            [ac_gen_scalar=no])

case ${ac_gen_scalar} in
  yes)
    AC_DEFINE([GENERIC_SCALAR],[1],[Use scalar data parallel loops])
  ;;
  *)
  ;;
esac

##################### Compiler dependent choices

#Strip any optional compiler arguments from nvcc call (eg -ccbin) for compiler comparison
CXXBASE=${CXX}
CXXTEST=${CXX}
if echo "${CXX}" | grep -q "nvcc"; then
  CXXTEST="nvcc"
fi

case ${CXXTEST} in
  nvcc)
#    CXX="nvcc -keep -v -x cu "
#    CXXLD="nvcc -v -link"
    CXX="${CXXBASE} -x cu "
    CXXLD="${CXXBASE} -link"
    CXXFLAGS="$CXXFLAGS -Xcompiler -fno-strict-aliasing --expt-extended-lambda --expt-relaxed-constexpr"
    if test $ac_openmp = yes; then
       CXXFLAGS="$CXXFLAGS -Xcompiler -fopenmp"
       LDFLAGS="$LDFLAGS -Xcompiler -fopenmp"
    fi
    ;;
  hipcc)
    CXXFLAGS="$CXXFLAGS -fno-strict-aliasing"
    CXXLD=${CXX}
    if test $ac_openmp = yes; then
       CXXFLAGS="$CXXFLAGS -fopenmp"
    fi
    ;;
  dpcpp)
    LDFLAGS="$LDFLAGS"
    CXXFLAGS="$CXXFLAGS"
    CXXLD=${CXX}
    ;;
  *)
    CXXLD=${CXX}
    CXXFLAGS="$CXXFLAGS -fno-strict-aliasing"
    ;;
esac

case ${ax_cv_cxx_compiler_vendor} in
  clang|gnu)
    case ${ac_SIMD} in
      GPU)
        AC_DEFINE([GPU_VEC],[1],[GPU vectorised])
        AC_DEFINE_UNQUOTED([GEN_SIMD_WIDTH],[$ac_gen_simd_width],
                           [generic SIMD vector width (in bytes)])
        SIMD_GEN_WIDTH_MSG=" (width= $ac_gen_simd_width)"
	SIMD_FLAGS='';;
      GPU-RRII)
        AC_DEFINE([GPU_RRII],[1],[GPU vectorised with RRRR / IIII layout])
        AC_DEFINE_UNQUOTED([GEN_SIMD_WIDTH],[$ac_gen_simd_width],
                           [generic SIMD vector width (in bytes)])
        SIMD_GEN_WIDTH_MSG=" (width= $ac_gen_simd_width)"
	SIMD_FLAGS='';;
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
        SIMD_FLAGS='-mavx512f -mavx512cd';;
      SKL)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics for SkyLake Xeon])
        SIMD_FLAGS='-march=skylake-avx512';;
      KNC)
        AC_DEFINE([IMCI],[1],[IMCI intrinsics for Knights Corner])
        SIMD_FLAGS='';;
      KNL)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics])
        AC_DEFINE([KNL],[1],[Knights landing processor])
        SIMD_FLAGS='-march=knl';;
      GEN)
        AC_DEFINE([GEN],[1],[generic vector code])
        AC_DEFINE_UNQUOTED([GEN_SIMD_WIDTH],[$ac_gen_simd_width],
                           [generic SIMD vector width (in bytes)])
        SIMD_GEN_WIDTH_MSG=" (width= $ac_gen_simd_width)"
        SIMD_FLAGS='';;
      A64FX)
        case ${ax_cv_cxx_compiler_vendor} in
          gnu)
            AC_DEFINE([A64FX],[1],[A64FX / 512-bit SVE VLS])
            SIMD_FLAGS='-march=armv8.2-a+sve -msve-vector-bits=512 -fno-gcse -DDSLASHINTRIN';;
          clang)
            AC_DEFINE([A64FX],[1],[A64FX / 512-bit SVE VLA])
            SIMD_FLAGS='-mcpu=a64fx -DARMCLANGCOMPAT -DDSLASHINTRIN';;
        esac;;
      NEONv8)
        AC_DEFINE([NEONV8],[1],[ARMv8 NEON])
        SIMD_FLAGS='-march=armv8-a';;
      QPX|BGQ)
        AC_DEFINE([QPX],[1],[QPX intrinsics for BG/Q])
        SIMD_FLAGS='';;
      *)
        AC_MSG_ERROR(["SIMD option ${ac_SIMD} not supported by the GCC/Clang compiler"]);;
    esac;;
  intel)
    case ${ac_SIMD} in
      GPU)
        AC_DEFINE([GPU_VEC],[1],[GPU vectorised ])
	SIMD_FLAGS='';;
      GPU-RRII)
        AC_DEFINE([GPU_RRII],[1],[GPU vectorised with RRRR / IIII layout])
        AC_DEFINE_UNQUOTED([GEN_SIMD_WIDTH],[$ac_gen_simd_width],
                           [generic SIMD vector width (in bytes)])
        SIMD_GEN_WIDTH_MSG=" (width= $ac_gen_simd_width)"
	SIMD_FLAGS='';;
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
        SIMD_FLAGS='-xcommon-avx512';;
      SKL)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics])
        SIMD_FLAGS='-xcore-avx512';;
      KNC)
        AC_DEFINE([IMCI],[1],[IMCI Intrinsics for Knights Corner])
        SIMD_FLAGS='';;
      KNL)
        AC_DEFINE([AVX512],[1],[AVX512 intrinsics for Knights Landing])
        AC_DEFINE([KNL],[1],[Knights landing processor])
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

###### PRECISION ALWAYS DOUBLE
AC_DEFINE([GRID_DEFAULT_PRECISION_DOUBLE],[1],[GRID_DEFAULT_PRECISION is DOUBLE] )

#########################################################
######################  GRID ALLOCATOR ALIGNMENT ##
#########################################################
AC_ARG_ENABLE([alloc-align],[AS_HELP_STRING([--enable-alloc-align=2MB|4k],[Alignment in bytes of GRID Allocator ])],[ac_ALLOC_ALIGN=${enable_alloc_align}],[ac_ALLOC_ALIGN=2MB])
case ${ac_ALLOC_ALIGN} in
    4k)
     AC_DEFINE([GRID_ALLOC_ALIGN],[(4096)],[GRID_ALLOC_ALIGN]);;
    2MB)
     AC_DEFINE([GRID_ALLOC_ALIGN],[(2*1024*1024)],[GRID_ALLOC_ALIGN]);;
    *);;
esac

AC_ARG_ENABLE([alloc-cache],[AS_HELP_STRING([--enable-alloc-cache ],[Cache a pool of recent "frees" to reuse])],[ac_ALLOC_CACHE=${enable_alloc_cache}],[ac_ALLOC_CACHE=yes])
case ${ac_ALLOC_CACHE} in
    yes)
     AC_DEFINE([ALLOCATION_CACHE],[1],[ALLOCATION_CACHE]);;
    *);;
esac


#########################################################
######################  set GPU device to rank in node ##
#########################################################
AC_ARG_ENABLE([setdevice],[AS_HELP_STRING([--enable-setdevice | --disable-setdevice],[Set GPU to rank in node with cudaSetDevice or similar])],[ac_SETDEVICE=${enable_SETDEVICE}],[ac_SETDEVICE=no])
case ${ac_SETDEVICE} in
    yes)
	echo ENABLE SET DEVICE
    ;;
    *)
     AC_DEFINE([GRID_DEFAULT_GPU],[1],[GRID_DEFAULT_GPU] )
   	echo DISABLE SET DEVICE
    ;;
esac

#########################################################
######################  Shared memory intranode #########
#########################################################
AC_ARG_ENABLE([shm],[AS_HELP_STRING([--enable-shm=shmopen|shmget|hugetlbfs|shmnone|nvlink|no|none],[Select SHM allocation technique])],[ac_SHM=${enable_shm}],[ac_SHM=no])

case ${ac_SHM} in

     shmopen)
     AC_DEFINE([GRID_MPI3_SHMOPEN],[1],[GRID_MPI3_SHMOPEN] )
     CXXFLAGS_CPY=$CXXFLAGS
     LDFLAGS_CPY=$LDFLAGS
     CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
     LDFLAGS="$AM_LDFLAGS $LDFLAGS"
     AC_SEARCH_LIBS([shm_unlink], [rt], [],[AC_MSG_ERROR("no library found for shm_unlink")])
     CXXFLAGS=$CXXFLAGS_CPY
     LDFLAGS=$LDFLAGS_CPY
     ;;

     shmget)
     AC_DEFINE([GRID_MPI3_SHMGET],[1],[GRID_MPI3_SHMGET] )
     ;;

     shmnone | no | none)
     AC_DEFINE([GRID_MPI3_SHM_NONE],[1],[GRID_MPI3_SHM_NONE] )
     ;;

     nvlink)
     AC_DEFINE([GRID_MPI3_SHM_NVLINK],[1],[GRID_MPI3_SHM_NVLINK] )
     ;;

     hugetlbfs)
     AC_DEFINE([GRID_MPI3_SHMMMAP],[1],[GRID_MPI3_SHMMMAP] )
     ;;

     *)
     AC_MSG_ERROR([${ac_SHM} unsupported --enable-shm option]);
     ;;
esac

######################  Shared base path for SHMMMAP
AC_ARG_ENABLE([shmpath],[AS_HELP_STRING([--enable-shmpath=path],[Select SHM mmap base path for hugetlbfs])],
	      [ac_SHMPATH=${enable_shmpath}],
	      [ac_SHMPATH=/var/lib/hugetlbfs/global/pagesize-2MB/])
AC_DEFINE_UNQUOTED([GRID_SHM_PATH],["$ac_SHMPATH"],[Path to a hugetlbfs filesystem for MMAPing])

############### force MPI in SMP
AC_ARG_ENABLE([shm-force-mpi],[AS_HELP_STRING([--enable-shm-force-mpi],[Force MPI within shared memory])],[ac_SHM_FORCE_MPI=${enable_shm_force_mpi}],[ac_SHM_FORCE_MPI=no])
case ${ac_SHM_FORCE_MPI} in
     yes)
        AC_DEFINE([GRID_SHM_FORCE_MPI],[1],[GRID_SHM_FORCE_MPI] )
      ;;
     *) ;;
esac
############### force MPI in SMP
AC_ARG_ENABLE([shm-fast-path],[AS_HELP_STRING([--enable-shm-fast-path],[Allow kernels to remote copy over intranode])],[ac_SHM_FAST_PATH=${enable_shm_fast_path}],[ac_SHM_FAST_PATH=no])
case ${ac_SHM_FAST_PATH} in
     yes)
        AC_DEFINE([SHM_FAST_PATH],[1],[SHM_FAST_PATH] )
      ;;
     *) ;;
esac

############### communication type selection
AC_ARG_ENABLE([comms-threads],[AS_HELP_STRING([--enable-comms-threads | --disable-comms-threads],[Use multiple threads in MPI calls])],[ac_COMMS_THREADS=${enable_comms_threads}],[ac_COMMS_THREADS=yes])

case ${ac_COMMS_THREADS} in
     yes)
        AC_DEFINE([GRID_COMMS_THREADING],[1],[GRID_COMMS_NONE] )
      ;;
     *) ;;
esac

############### communication type selection
AC_ARG_ENABLE([comms],[AS_HELP_STRING([--enable-comms=none|mpi|mpi-auto],[Select communications])],[ac_COMMS=${enable_comms}],[ac_COMMS=none])


case ${ac_COMMS} in
     none)
        AC_DEFINE([GRID_COMMS_NONE],[1],[GRID_COMMS_NONE] )
        comms_type='none'
     ;;
     mpi*)
        AC_DEFINE([GRID_COMMS_MPI3],[1],[GRID_COMMS_MPI3] )
        comms_type='mpi3'
     ;;
     *)
        AC_MSG_ERROR([${ac_COMMS} unsupported --enable-comms option]);
     ;;
esac
case ${ac_COMMS} in
    *-auto)
        LX_FIND_MPI
	echo MPI_CXXFLAGS $MPI_CXXFLAGS
	echo MPI_CXXLDFLAGS $MPI_CXXLDFLAGS
	echo MPI_CFLAGS $MPI_CFLAGS
        AM_CXXFLAGS="$MPI_CXXFLAGS $AM_CXXFLAGS"
        AM_CFLAGS="$MPI_CFLAGS $AM_CFLAGS"
        AM_LDFLAGS="`echo $MPI_CXXLDFLAGS | sed -E 's/-l@<:@^ @:>@+//g'` $AM_LDFLAGS"
        LIBS="`echo $MPI_CXXLDFLAGS | sed -E 's/-L@<:@^ @:>@+//g'` $LIBS";;
    *)
        ;;
esac

AM_CONDITIONAL(BUILD_COMMS_MPI3,  [ test "${comms_type}X" == "mpi3X" ] )
AM_CONDITIONAL(BUILD_COMMS_NONE,  [ test "${comms_type}X" == "noneX" ])

############### RNG selection
AC_ARG_ENABLE([rng],[AS_HELP_STRING([--enable-rng=ranlux48|mt19937|sitmo],[\
	            Select Random Number Generator to be used])],\
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
AC_ARG_ENABLE([timers],[AS_HELP_STRING([--enable-timers],[\
	            Enable system dependent high res timers])],\
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
AC_ARG_ENABLE([chroma],[AS_HELP_STRING([--enable-chroma],[Expect chroma compiled under c++14 ])],ac_CHROMA=yes,ac_CHROMA=no)

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
GRID_CXX="$CXX"
GRID_CXXLD="$CXXLD"
GRID_CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
GRID_LDFLAGS="$AM_LDFLAGS $LDFLAGS"
GRID_LIBS=$LIBS
GRID_SHORT_SHA=`git rev-parse --short HEAD`
GRID_SHA=`git rev-parse HEAD`
GRID_BRANCH=`git rev-parse --abbrev-ref HEAD`
AM_CXXFLAGS="-I${abs_srcdir} $AM_CXXFLAGS"
AM_CFLAGS="-I${abs_srcdir} $AM_CFLAGS"
AM_LDFLAGS="-L${cwd}/Grid $AM_LDFLAGS"
AC_SUBST([CXXLD])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([GRID_CXX])
AC_SUBST([GRID_CXXLD])
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
Nc                          : ${ac_Nc}
SIMD                        : ${ac_SIMD}${SIMD_GEN_WIDTH_MSG}
Threading                   : ${ac_openmp}
Acceleration                : ${ac_ACCELERATOR}
Unified virtual memory      : ${ac_UNIFIED}
Communications type         : ${comms_type}
Shared memory allocator     : ${ac_SHM}
Shared memory mmap path     : ${ac_SHMPATH}
Default precision           : ${ac_PRECISION}
Software FP16 conversion    : ${ac_SFW_FP16}
RNG choice                  : ${ac_RNG}
GMP                         : `if test "x$have_gmp" = xtrue; then echo yes; else echo no; fi`
LAPACK                      : ${ac_LAPACK}
FFTW                        : `if test "x$have_fftw" = xtrue; then echo yes; else echo no; fi`
LIME (ILDG support)         : `if test "x$have_lime" = xtrue; then echo yes; else echo no; fi`
HDF5                        : `if test "x$have_hdf5" = xtrue; then echo yes; else echo no; fi`
build DOXYGEN documentation : `if test "$DX_FLAG_doc" = '1'; then echo yes; else echo no; fi`
Sp2n                        : ${ac_ENABLE_SP}
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
AC_CONFIG_FILES(Grid/Makefile)
AC_CONFIG_FILES(HMC/Makefile)
AC_CONFIG_FILES(tests/Makefile)
AC_CONFIG_FILES(tests/IO/Makefile)
AC_CONFIG_FILES(tests/core/Makefile)
AC_CONFIG_FILES(tests/debug/Makefile)
AC_CONFIG_FILES(tests/forces/Makefile)
AC_CONFIG_FILES(tests/hmc/Makefile)
AC_CONFIG_FILES(tests/solver/Makefile)
AC_CONFIG_FILES(tests/lanczos/Makefile)
AC_CONFIG_FILES(tests/smearing/Makefile)
AC_CONFIG_FILES(tests/qdpxx/Makefile)
AC_CONFIG_FILES(tests/testu01/Makefile)
AC_CONFIG_FILES(tests/sp2n/Makefile)
AC_CONFIG_FILES(benchmarks/Makefile)
AC_CONFIG_FILES(examples/Makefile)
AC_OUTPUT

echo ""
cat grid.configure.summary
echo ""
