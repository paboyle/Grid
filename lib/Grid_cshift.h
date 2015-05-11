#ifndef _GRID_CSHIFT_H_
#define _GRID_CSHIFT_H_

#include <cshift/Grid_cshift_common.h>

#ifdef GRID_COMMS_NONE
#include <cshift/Grid_cshift_none.h>
#endif

#ifdef GRID_COMMS_MPI
#include <cshift/Grid_cshift_mpi.h>
#endif 
#endif
