#!/bin/bash

#export NUMA_MAP=(2 2 2 3 3 3 2 2 2 3 3 3 )
#export NUMA_MAP=(0 0 1 1 0 0 1 1 0 0 1 1);
#export  GPU_MAP=(0.0 0.1 3.0 3.1 1.0 1.1 4.0 4.1 2.0 2.1 5.0 5.1)

export NUMA_MAP=(0 0 0 0 0 0 1 1 1 1 1 1 );
export  GPU_MAP=(0.0 1.0 2.0 3.0 4.0 5.0 0.1 1.1 2.1 3.1 4.1 5.1 )

export NUMA=${NUMA_MAP[$PALS_LOCAL_RANKID]}
export gpu_id=${GPU_MAP[$PALS_LOCAL_RANKID]}
  
unset EnableWalkerPartition
export EnableImplicitScaling=0
export ZE_AFFINITY_MASK=$gpu_id
export ONEAPI_DEVICE_FILTER=gpu,level_zero

export SYCL_PI_LEVEL_ZERO_DEVICE_SCOPE_EVENTS=0
export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=1
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE=0:5
#export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE=0:2
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE_FOR_D2D_COPY=1
#export SYCL_PI_LEVEL_ZERO_USM_RESIDENT=1

echo "rank $PALS_RANKID ; local rank $PALS_LOCAL_RANKID ; ZE_AFFINITY_MASK=$ZE_AFFINITY_MASK ; NUMA $NUMA "

if [ $PALS_RANKID = "0" ]
then
#    numactl -m $NUMA -N $NUMA onetrace --chrome-device-timeline  "$@"
#    numactl -m $NUMA -N $NUMA unitrace --chrome-kernel-logging --chrome-mpi-logging --chrome-sycl-logging --demangle "$@"
    numactl -m $NUMA -N $NUMA  "$@"
else 
    numactl -m $NUMA -N $NUMA  "$@"
fi
