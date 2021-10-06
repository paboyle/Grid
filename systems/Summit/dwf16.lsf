#!/bin/bash
#BSUB -P LGT104
#BSUB -W 2:00
#BSUB -nnodes 16
#BSUB -J DWF

export OMP_NUM_THREADS=6
export PAMI_IBV_ADAPTER_AFFINITY=1
export PAMI_ENABLE_STRIPING=1
export OPT="--comms-concurrent --comms-overlap "

APP="./benchmarks/Benchmark_comms_host_device  --mpi 4.4.4.3 "
jsrun --nrs 16 -a6 -g6 -c42 -dpacked -b packed:7 --latency_priority gpu-cpu --smpiargs=-gpu $APP > comms.16node.log

APP="./benchmarks/Benchmark_dwf_fp32 --grid 96.96.96.72 --mpi 4.4.4.3 --shm 2048 --shm-force-mpi 1 --device-mem 8000 --shm-force-mpi 1 $OPT "
jsrun --nrs 16 -a6 -g6 -c42 -dpacked -b packed:7 --latency_priority gpu-cpu --smpiargs=-gpu $APP > dwf.16node.24.log

APP="./benchmarks/Benchmark_dwf_fp32 --grid 128.128.128.96 --mpi 4.4.4.3 --shm 2048 --shm-force-mpi 1 --device-mem 8000 --shm-force-mpi 1 $OPT "
jsrun --nrs 16 -a6 -g6 -c42 -dpacked -b packed:7 --latency_priority gpu-cpu --smpiargs=-gpu $APP > dwf.16node.32.log






