#ifndef _GRID_CSHIFT_H_
#define _GRID_CSHIFT_H_

#include <cshift/Cshift_common.h>

#ifdef GRID_COMMS_NONE
#include <cshift/Cshift_none.h>
#endif

#ifdef GRID_COMMS_MPI
#include <cshift/Cshift_mpi.h>
#endif 
#endif
