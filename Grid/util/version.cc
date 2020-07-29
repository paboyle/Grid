#include <iostream>
#include "Version.h"
namespace Grid {
  void printHash(){
#ifdef GITHASH
    std::cout << "Current Grid git commit hash=" << GITHASH << std::endl;
#else
    std::cout << "Current Grid git commit hash is undefined. Check makefile." << std::endl;
#endif
#undef GITHASH
}
}
