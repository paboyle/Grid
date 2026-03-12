#pragma once

#ifndef BUILD_FERMION_INSTANTIATIONS
#include <iostream>

int main(void) {
  std::cout << "This build of Grid was configured to exclude fermion instantiations, "
	    << "which this benchmark relies on. "
	    << "Please reconfigure and rebuild Grid with --enable-fermion-instantiations"
	    << "to run this benchmark."
	    << std::endl;
  return 1;
}
#endif
