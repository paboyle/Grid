#!/bin/bash

#export NUMA_MAP=(2 2 2 3 3 3 2 2 2 3 3 3 )
#export NUMA_MAP=(0 0 1 1 0 0 1 1 0 0 1 1);
#export  GPU_MAP=(0.0 0.1 3.0 3.1 1.0 1.1 4.0 4.1 2.0 2.1 5.0 5.1)

export NUMA_PMAP=(0 0 0 1 1 1 0 0 0 1 1 1 );
export NUMA_HMAP=(2 2 2 3 3 3 3 2 2 2 2 3 3 3 );
export  GPU_MAP=(0.0 1.0 2.0 3.0 4.0 5.0 0.1 1.1 2.1 3.1 4.1 5.1 )

export NUMAP=${NUMA_PMAP[$PALS_LOCAL_RANKID]}
export NUMAH=${NUMA_HMAP[$PALS_LOCAL_RANKID]}
export gpu_id=${GPU_MAP[$PALS_LOCAL_RANKID]}
  
unset EnableWalkerPartition
export EnableImplicitScaling=0
export ZE_AFFINITY_MASK=$gpu_id
export ONEAPI_DEVICE_FILTER=gpu,level_zero

export SYCL_PI_LEVEL_ZERO_DEVICE_SCOPE_EVENTS=0
export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=1
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE=0:4
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE_FOR_D2D_COPY=1
#export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE=0:2
#export SYCL_PI_LEVEL_ZERO_USM_RESIDENT=1

#export MPI_BUF_NUMA=$NUMAH

echo "rank $PALS_RANKID ; local rank $PALS_LOCAL_RANKID ; ZE_AFFINITY_MASK=$ZE_AFFINITY_MASK ; NUMA $NUMA "

if [ $PALS_RANKID = "0" ]
then
#    numactl -p $NUMAP -N $NUMAP unitrace --chrome-kernel-logging --chrome-mpi-logging --chrome-sycl-logging --demangle "$@"
    numactl -p $NUMAP -N $NUMAP  "$@"
else 
    numactl -p $NUMAP -N $NUMAP  "$@"
fi
