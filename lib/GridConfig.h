/* lib/GridConfig.h.  Generated from GridConfig.h.in by configure.  */
/* lib/GridConfig.h.in.  Generated from configure.ac by autoheader.  */

/* AVX Intrinsics */
/* #undef AVX1 */

/* AVX2 Intrinsics */
/* #undef AVX2 */

/* AVX512 Intrinsics for Knights Corner */
/* #undef AVX512 */

/* EMPTY_SIMD only for DEBUGGING */
/* #undef EMPTY_SIMD */

/* GRID_COMMS_MPI */
/* #undef GRID_COMMS_MPI */

/* GRID_COMMS_NONE */
#define GRID_COMMS_NONE 1

/* GRID_DEFAULT_PRECISION is DOUBLE */
/* #undef GRID_DEFAULT_PRECISION_DOUBLE */

/* GRID_DEFAULT_PRECISION is SINGLE */
#define GRID_DEFAULT_PRECISION_SINGLE 1

/* Support Altivec instructions */
/* #undef HAVE_ALTIVEC */

/* Support AVX (Advanced Vector Extensions) instructions */
/* #undef HAVE_AVX */

/* Support AVX2 (Advanced Vector Extensions 2) instructions */
/* #undef HAVE_AVX2 */

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

/* Support FMA3 (Fused Multiply-Add) instructions */
/* #undef HAVE_FMA */

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <gmp.h> header file. */
#define HAVE_GMP_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `gmp' library (-lgmp). */
#define HAVE_LIBGMP 1

/* Define to 1 if you have the `mpfr' library (-lmpfr). */
#define HAVE_LIBMPFR 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <malloc/malloc.h> header file. */
/* #undef HAVE_MALLOC_MALLOC_H */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Support mmx instructions */
#define HAVE_MMX /**/

/* Define to 1 if you have the <mm_malloc.h> header file. */
#define HAVE_MM_MALLOC_H 1

/* Support SSE (Streaming SIMD Extensions) instructions */
#define HAVE_SSE /**/

/* Support SSE2 (Streaming SIMD Extensions 2) instructions */
#define HAVE_SSE2 /**/

/* Support SSE3 (Streaming SIMD Extensions 3) instructions */
#define HAVE_SSE3 /**/

/* Support SSSE4.1 (Streaming SIMD Extensions 4.1) instructions */
#define HAVE_SSE4_1 /**/

/* Support SSSE4.2 (Streaming SIMD Extensions 4.2) instructions */
#define HAVE_SSE4_2 /**/

/* Support SSSE3 (Supplemental Streaming SIMD Extensions 3) instructions */
#define HAVE_SSSE3 /**/

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

/* NEON ARMv8 Experimental support */
/* #undef NEONv8 */

/* Name of package */
#define PACKAGE "grid"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "paboyle@ph.ed.ac.uk"

/* Define to the full name of this package. */
#define PACKAGE_NAME "Grid"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "Grid 1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "grid"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0"

/* SSE4 Intrinsics */
#define SSE4 1

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "1.0"

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
