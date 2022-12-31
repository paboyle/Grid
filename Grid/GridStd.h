#ifndef GRID_STD_H
#define GRID_STD_H

///////////////////
// Std C++ dependencies
///////////////////
#include <cassert>
#include <complex>
#include <memory>
#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <signal.h>
#include <ctime>
#include <sys/time.h>
#include <chrono>
#include <zlib.h>

///////////////////
// Grid config
///////////////////
#include "Config.h"

#ifdef TOFU
#undef GRID_COMMS_THREADS
#endif
#endif /* GRID_STD_H */
