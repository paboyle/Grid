#pragma once 


#define COMPILER_FATAL     
#define COMPILER_GCC_WARN  "Compiler" COMPILER_STRING " is known to not work with Grid under -O3 due to compiler bugs"

#if defined(__clang__)

  #if __clang_major < 3
    #error "This clang++ version is known to not work with Grid due to compiler bugs"
  #endif

  #if __clang_major == 3
    #if __clang_minor < 5
    #error "This clang++ version is known to not work with Grid due to compiler bugs"
    #endif
  #endif 

#endif

// Intel compiler *ALSO* has __GNUC__ defined so must if/else GCC checks
#ifdef __INTEL_COMPILER

#if __INTEL_COMPILER < 1600
  #error "This icpc version is known to not work with Grid due to compiler bugs"
#endif

#else 

  #if defined(__GNUC__)

    #if __GNUC__ < 4 
      #error "g++ prior to version 4 is known to not work with Grid due to compiler bugs"
    #endif

    #if __GNUC__ == 4 
      #if __GNUC_MINOR__ != 9 
      #error "g++ 4.9 is the only gcc-4.x version known to work with Grid due to compiler bugs"
      #endif
    #endif

    #if __GNUC__ == 5
      #warning "g++ version 5 is known to not work with Grid due to compiler bugs under -O3 : ensure you run make check"
    #endif

    #if __GNUC__ == 6
      #if __GNUC_MINOR__ < 3
        #warning "This g++6.3 is the first recent g++ version known to work with Grid: ensure you run make check"
      #endif
    #endif

  #endif

#endif

