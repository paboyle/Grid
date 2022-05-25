#!/bin/bash
#SBATCH -A m3886_g
#SBATCH -C gpu
#SBATCH -q debug
#SBATCH -t 0:20:00
#SBATCH -c 32
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --exclusive
#SBATCH --gpu-bind=none

export SLURM_CPU_BIND="cores"
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_RDMA_ENABLED_CUDA=1
export MPICH_GPU_IPC_ENABLED=1
export MPICH_GPU_EAGER_REGISTER_HOST_MEM=0
export MPICH_GPU_NO_ASYNC_MEMCPY=0
#export MPICH_SMP_SINGLE_COPY_MODE=CMA

OPT="--comms-overlap --shm-mpi 1"
VOL=64.64.32.32
srun ./benchmarks/Benchmark_dwf_fp32 --mpi 2.2.1.1 --grid $VOL --accelerator-threads 8 --shm 2048 $OPT
#srun ./benchmarks/Benchmark_dwf_fp32 --mpi 2.1.1.4 --grid $VOL --accelerator-threads 8 --shm 2048 $OPT
#srun ./benchmarks/Benchmark_dwf_fp32 --mpi 1.1.1.8 --grid $VOL --accelerator-threads 8 --shm 2048 $OPT

