/* lib/Config.h.  Generated from Config.h.in by configure.  */
/* lib/Config.h.in.  Generated from configure.ac by autoheader.  */

/* AVX intrinsics */
#define AVX1 1

/* AVX2 intrinsics */
/* #undef AVX2 */

/* AVX512 intrinsics for Knights Landing */
/* #undef AVX512 */

/* AVX intrinsics with FMA3 */
/* #undef AVXFMA */

/* AVX intrinsics with FMA4 */
/* #undef AVXFMA4 */

/* vendor of C++ compiler that will compile the code */
#define CXX_COMP_VENDOR "gnu"

/* generic vector code */
/* #undef GEN */

/* generic vector code */
/* #undef GENERIC_VEC */

/* generic SIMD vector width (in bytes) */
/* #undef GEN_SIMD_WIDTH */

/* GRID_COMMS_MPI3 */
#define GRID_COMMS_MPI3 1

/* GRID_COMMS_NONE */
/* #undef GRID_COMMS_NONE */

/* GRID_DEFAULT_PRECISION is DOUBLE */
#define GRID_DEFAULT_PRECISION_DOUBLE 1

/* GRID_DEFAULT_PRECISION is SINGLE */
/* #undef GRID_DEFAULT_PRECISION_SINGLE */

/* GRID_MPI3_SHMGET */
/* #undef GRID_MPI3_SHMGET */

/* GRID_MPI3_SHMMMAP */
/* #undef GRID_MPI3_SHMMMAP */

/* GRID_MPI3_SHMOPEN */
#define GRID_MPI3_SHMOPEN 1

/* GRID_MPI3_SHM_NONE */
/* #undef GRID_MPI3_SHM_NONE */

/* First touch numa locality */
/* #undef GRID_NUMA */

/* Path to a hugetlbfs filesystem for MMAPing */
#define GRID_SHM_PATH "/var/lib/hugetlbfs/global/pagesize-2MB/"

/* version of g++ that will compile the code */
#define GXX_VERSION "6.4.0"

/* define if the compiler supports basic C++11 syntax */
/* #undef HAVE_CXX11 */

/* Define to 1 if you have the declaration of `be64toh', and to 0 if you
   don't. */
#define HAVE_DECL_BE64TOH 1

/* Define to 1 if you have the declaration of `ntohll', and to 0 if you don't.
   */
#define HAVE_DECL_NTOHLL 0

/* Define to 1 if you have the <endian.h> header file. */
#define HAVE_ENDIAN_H 1

/* Define to 1 if you have the <execinfo.h> header file. */
#define HAVE_EXECINFO_H 1

/* Define to 1 if you have the `FFTW' library */
/* #undef HAVE_FFTW */

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the `HDF5' library */
/* #undef HAVE_HDF5 */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `GMP' library */
/* #undef HAVE_LIBGMP */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `MPFR' library */
/* #undef HAVE_LIBMPFR */

/* Define to 1 if you have the `LIBNUMA' library */
#define HAVE_LIBNUMA 1

/* Define to 1 if you have the `stdc++' library (-lstdc++). */
#define HAVE_LIBSTDC__ 1

/* Define to 1 if you have the `LIME' library */
/* #undef HAVE_LIME */

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <malloc/malloc.h> header file. */
/* #undef HAVE_MALLOC_MALLOC_H */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <mm_malloc.h> header file. */
#define HAVE_MM_MALLOC_H 1

/* Define to 1 if you have MPI libs and headers. */
#define HAVE_MPI 1

/* Define to 1 if you have the <numaif.h> header file. */
#define HAVE_NUMAIF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `LIBZ' library */
#define HAVE_ZLIB 1

/* IMCI Intrinsics for Knights Corner */
/* #undef IMCI */

/* Knights landing processor */
/* #undef KNL */

/* ARMv8 NEON */
/* #undef NEONV8 */

/* Name of package */
#define PACKAGE "Grid"

/* Define to the address where bug reports for this package should be sent. */
#define GRID_BUGREPORT "https://github.com/paboyle/Grid"

/* Define to the full name of this package. */
#define GRID_NAME "Grid"

/* Define to the full name and version of this package. */
#define GRID_STRING "Grid 0.7.0"

/* Define to the one symbol short name of this package. */
#define GRID_TARNAME "Grid"

/* Define to the home page for this package. */
#define GRID_URL ""

/* Define to the version of this package. */
#define GRID_VERSION "0.7.0"

/* QPX intrinsics for BG/Q */
/* #undef QPX */

/* RNG_MT19937 */
/* #undef RNG_MT19937 */

/* RNG_RANLUX */
/* #undef RNG_RANLUX */

/* RNG_SITMO */
#define RNG_SITMO 1

/* software conversion to fp16 */
#define SFW_FP16 1

/* SSE4 intrinsics */
/* #undef SSE4 */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* compile ZMM test */
/* #undef TEST_ZMM */

/* TIMERS_OFF */
/* #undef TIMERS_OFF */

/* TIMERS_ON */
#define TIMERS_ON 1

/* use LAPACK */
/* #undef USE_LAPACK */

/* Define to 1 if you use the Intel MKL */
/* #undef USE_MKL */

/* Version number of package */
#define VERSION "0.7.0"

/* Define for Solaris 2.5.1 so the uint32_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT32_T */

/* Define for Solaris 2.5.1 so the uint64_t typedef from <sys/synch.h>,
   <pthread.h>, or <semaphore.h> is not used. If the typedef were allowed, the
   #define below would cause a syntax error. */
/* #undef _UINT64_T */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to the type of an unsigned integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint32_t */

/* Define to the type of an unsigned integer type of width exactly 64 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint64_t */
