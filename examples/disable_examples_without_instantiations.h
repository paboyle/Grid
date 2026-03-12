#pragma once

#ifndef BUILD_FERMION_INSTANTIATIONS
#include <iostream>

int main(void) {
  std::cout << "This build of Grid was configured to exclude fermion instantiations, "
	    << "which this example relies on. "
	    << "Please reconfigure and rebuild Grid with --enable-fermion-instantiations"
	    << "to run this example."
	    << std::endl;
  return 1;
}
#endif
