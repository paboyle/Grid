dnl @synopsis AX_GCC_X86_CPUID(OP)
dnl
dnl @summary run x86 cpuid instruction OP using gcc inline assembler
dnl
dnl On Pentium and later x86 processors, with gcc or a compiler that
dnl has a compatible syntax for inline assembly instructions, run a
dnl small program that executes the cpuid instruction with input OP.
dnl This can be used to detect the CPU type.
dnl
dnl On output, the values of the eax, ebx, ecx, and edx registers are
dnl stored as hexadecimal strings as "eax:ebx:ecx:edx" in the cache
dnl variable ax_cv_gcc_x86_cpuid_OP.
dnl
dnl If the cpuid instruction fails (because you are running a
dnl cross-compiler, or because you are not using gcc, or because you
dnl are on a processor that doesn't have this instruction),
dnl ax_cv_gcc_x86_cpuid_OP is set to the string "unknown".
dnl
dnl This macro mainly exists to be used in AX_GCC_ARCHFLAG.
dnl
dnl @category Misc
dnl @author Steven G. Johnson <stevenj@alum.mit.edu> and Matteo Frigo.
dnl @version 2005-05-30
dnl @license GPLWithACException

AC_DEFUN([AX_GCC_X86_CPUID],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
AC_CACHE_CHECK(for x86 cpuid $1 output, ax_cv_gcc_x86_cpuid_$1,
 [AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <stdio.h>], [
     int op = $1, eax, ebx, ecx, edx;
     FILE *f;
      __asm__("cpuid"
        : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx)
        : "a" (op));
     f = fopen("conftest_cpuid", "w"); if (!f) return 1;
     fprintf(f, "%x:%x:%x:%x\n", eax, ebx, ecx, edx);
     fclose(f);
     return 0;
])],
     [ax_cv_gcc_x86_cpuid_$1=`cat conftest_cpuid`; rm -f conftest_cpuid],
     [ax_cv_gcc_x86_cpuid_$1=unknown; rm -f conftest_cpuid],
     [ax_cv_gcc_x86_cpuid_$1=unknown])])
AC_LANG_POP([C])
])
