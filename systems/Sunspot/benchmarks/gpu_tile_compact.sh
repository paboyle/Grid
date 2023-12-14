#!/bin/bash

display_help() {
  echo " Will map gpu tile to rank in compact and then round-robin fashion"
  echo " Usage (only work for one node of ATS/PVC):"
  echo "   mpiexec --np N gpu_tile_compact.sh ./a.out"
  echo
  echo " Example 3 GPU of 2 Tiles with 7 Ranks:"
  echo "   0 Rank 0.0"
  echo "   1 Rank 0.1"
  echo "   2 Rank 1.0"
  echo "   3 Rank 1.1"
  echo "   4 Rank 2.0"
  echo "   5 Rank 2.1"
  echo "   6 Rank 0.0"
  echo
  echo " Hacked together by apl@anl.gov, please contact if bug found"
  exit 1
}

#This give the exact GPU count i915 knows about and I use udev to only enumerate the devices with physical presence.
#works? num_gpu=$(/usr/bin/udevadm info /sys/module/i915/drivers/pci\:i915/* |& grep -v Unknown | grep -c "P: /devices")
num_gpu=6
num_tile=2

if [ "$#" -eq 0 ] || [ "$1" == "--help" ] || [ "$1" == "-h" ] || [ "$num_gpu" = 0 ]; then
  display_help
fi

gpu_id=$(( (PALS_LOCAL_RANKID / num_tile ) % num_gpu ))
tile_id=$((PALS_LOCAL_RANKID % num_tile))

unset EnableWalkerPartition
export EnableImplicitScaling=0
export ZE_ENABLE_PCI_ID_DEVICE_ORDER=1
export ZE_AFFINITY_MASK=$gpu_id.$tile_id
export ONEAPI_DEVICE_FILTER=gpu,level_zero
export SYCL_PI_LEVEL_ZERO_DEVICE_SCOPE_EVENTS=0
export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=1
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE=0:2
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE_FOR_D2D_COPY=1
#export SYCL_PI_LEVEL_ZERO_USM_RESIDENT=1

echo "rank $PALS_RANKID ; local rank $PALS_LOCAL_RANKID ; ZE_AFFINITY_MASK=$ZE_AFFINITY_MASK"

if [ $PALS_LOCAL_RANKID = 0 ]
then
    onetrace --chrome-device-timeline "$@"
#    "$@"
else
"$@"
fi
