#ifndef _GRID_CSHIFT_H_
#define _GRID_CSHIFT_H_
#include <Grid_cshift_common.h>

#ifdef GRID_COMMS_NONE
#include <Grid_cshift_none.h>
#endif

#ifdef GRID_COMMS_FAKE
#include <Grid_cshift_fake.h>
#endif

#ifdef GRID_COMMS_MPI
#include <Grid_cshift_mpi.h>
#endif 
#endif
