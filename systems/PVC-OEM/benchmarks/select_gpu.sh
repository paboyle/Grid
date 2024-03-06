#!/bin/bash

num_tile=2
gpu_id=$(( (MPI_LOCALRANKID / num_tile ) ))
tile_id=$((MPI_LOCALRANKID % num_tile))

export ZE_AFFINITY_MASK=$gpu_id.$tile_id

echo "local rank $MPI_LOCALRANKID ; ZE_AFFINITY_MASK=$ZE_AFFINITY_MASK"

"$@"

