#!/bin/bash

export NUMA_MAP=(2 2 2 3 3 3 2 2 2 3 3 3 )
#export NUMA_MAP=(0 0 0 1 1 1 0 0 0 1 1 1 )
export NUMA_PMAP=(0 0 0 1 1 1 0 0 0 1 1 1 )
export  NIC_MAP=(0 1 2 4 5 6 0 1 2 4 5 6 )
export  GPU_MAP=(0 1 2 3 4 5 0 1 2 3 4 5 )
export TILE_MAP=(0 0 0 0 0 0 1 1 1 1 1 1 )

export NUMA=${NUMA_MAP[$PALS_LOCAL_RANKID]}
export NUMAP=${NUMA_PMAP[$PALS_LOCAL_RANKID]}
export NIC=${NIC_MAP[$PALS_LOCAL_RANKID]}
export gpu_id=${GPU_MAP[$PALS_LOCAL_RANKID]}
export tile_id=${TILE_MAP[$PALS_LOCAL_RANKID]}
  
#export GRID_MPICH_NIC_BIND=$NIC
#export ONEAPI_DEVICE_SELECTOR=level_zero:$gpu_id.$tile_id

unset EnableWalkerPartition
export EnableImplicitScaling=0
export ZE_AFFINITY_MASK=$gpu_id.$tile_id
export ONEAPI_DEVICE_FILTER=gpu,level_zero

#export ZE_ENABLE_PCI_ID_DEVICE_ORDER=1
#export SYCL_PI_LEVEL_ZERO_DEVICE_SCOPE_EVENTS=0
#export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=1
#export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE=0:2
#export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE_FOR_D2D_COPY=1
#export SYCL_PI_LEVEL_ZERO_USM_RESIDENT=1

#echo "rank $PALS_RANKID ; local rank $PALS_LOCAL_RANKID ; ZE_AFFINITY_MASK=$ZE_AFFINITY_MASK ; NUMA $NUMA "

numactl -m $NUMA -N $NUMAP  "$@"
