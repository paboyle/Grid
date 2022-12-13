#!/bin/sh

export ZE_AFFINITY_MASK=0.$MPI_LOCALRANKID

echo Ranke $MPI_LOCALRANKID ZE_AFFINITY_MASK is $ZE_AFFINITY_MASK


if [ $MPI_LOCALRANKID = "0" ] 
then
#  ~psteinbr/build_pti/ze_tracer -h $@
  onetrace --chrome-device-timeline $@
else
  $@
fi
