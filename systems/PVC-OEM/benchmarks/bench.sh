#!/bin/bash

export EnableImplicitScaling=0
export ZE_ENABLE_PCI_ID_DEVICE_ORDER=1
export ZE_AFFINITY_MASK=$gpu_id.$tile_id
export ONEAPI_DEVICE_FILTER=gpu,level_zero
export SYCL_PI_LEVEL_ZERO_DEVICE_SCOPE_EVENTS=0
export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=1
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE=0:2
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE_FOR_D2D_COPY=1

mpiexec -launcher ssh -n 1 -host localhost ./select_gpu.sh ./Benchmark_dwf_fp32 --mpi 1.1.1.1 --grid 32.32.32.32 --accelerator-threads 16 --shm-mpi 1 --shm 2048 --device-mem 32768 | tee 1tile.log
mpiexec -launcher ssh -n 2 -host localhost ./select_gpu.sh ./Benchmark_dwf_fp32 --mpi 1.1.1.2 --grid 32.32.32.64 --accelerator-threads 16 --shm-mpi 1 --shm 2048 --device-mem 32768 | tee 2tile.log

#mpiexec -launcher ssh -n 4 -host localhost ./select_gpu.sh ./Benchmark_dwf_fp32 --mpi 1.1.2.2 --grid 16.16.64.64 --accelerator-threads 16 --shm-mpi 0 --shm 2048 --device-mem 32768 | tee 4tile.log
#mpiexec -launcher ssh -n 8 -host localhost ./select_gpu.sh ./Benchmark_dwf_fp32 --mpi 1.1.2.4 --grid 16.16.64.128 --accelerator-threads 16 --shm-mpi 0 --shm 2048 --device-mem 32768 | tee 8tile.log


