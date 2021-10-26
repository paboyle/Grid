#!/bin/bash
#BSUB -P LGT104
#BSUB -W 2:00
#BSUB -nnodes 4
#BSUB -J DWF

export OMP_NUM_THREADS=6
export PAMI_IBV_ADAPTER_AFFINITY=1
export PAMI_ENABLE_STRIPING=1
export OPT="--comms-concurrent --comms-overlap "
#export GRID_ALLOC_NCACHE_LARGE=1
export APP="./benchmarks/Benchmark_comms_host_device  --mpi 2.2.2.3 "
jsrun --nrs 4 -a6 -g6 -c42 -dpacked -b packed:7 --latency_priority gpu-cpu --smpiargs=-gpu $APP > comms.4node

APP="./benchmarks/Benchmark_dwf_fp32 --grid 48.48.48.72 --mpi 2.2.2.3 --shm 2048 --shm-force-mpi 1 --device-mem 8000 --shm-force-mpi 1 $OPT "
jsrun --nrs 4 -a6 -g6 -c42 -dpacked -b packed:7 --latency_priority gpu-cpu --smpiargs=-gpu $APP > dwf.24.4node

APP="./benchmarks/Benchmark_dwf_fp32 --grid 64.64.64.96 --mpi 2.2.2.3 --shm 2048 --shm-force-mpi 1 --device-mem 8000 --shm-force-mpi 1 $OPT "
jsrun --nrs 4 -a6 -g6 -c42 -dpacked -b packed:7 --latency_priority gpu-cpu --smpiargs=-gpu $APP > dwf.32.4node






