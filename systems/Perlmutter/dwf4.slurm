#!/bin/bash
#SBATCH -A m3886_g
#SBATCH -C gpu
#SBATCH -q debug
#SBATCH -t 0:10:00
#SBATCH -n 4
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH --exclusive
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=none

export SLURM_CPU_BIND="cores"
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_RDMA_ENABLED_CUDA=0
export MPICH_GPU_IPC_ENABLED=0
export MPICH_GPU_EAGER_REGISTER_HOST_MEM=0
export MPICH_GPU_NO_ASYNC_MEMCPY=1
OPT="--comms-overlap --comms-concurrent --shm-mpi 1"
srun ./benchmarks/Benchmark_ITT --mpi 2.1.1.2 --grid 64.64.64.64 --accelerator-threads 8 --shm 2048 $OPT > ITT.log
