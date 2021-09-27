#!/bin/bash
#SBATCH -A mp13
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 0:20:00
#SBATCH -n 16
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH --exclusive
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=map_gpu:0,1,2,3

export SLURM_CPU_BIND="cores"
export MPICH_RDMA_ENABLED_CUDA=1
export MPICH_GPU_SUPPORT_ENABLED=1
srun ./benchmarks/Benchmark_comms_host_device --mpi 2.2.2.2  --accelerator-threads 8 > comms.4node

OPT="--comms-overlap --comms-concurrent --shm-mpi 0"
srun ./benchmarks/Benchmark_dwf_fp32 --mpi 2.2.2.2 --grid 64.64.64.64 --accelerator-threads 8 --shm 2048 $OPT > dwf.64.64.64.64.4node.opt0
srun ./benchmarks/Benchmark_dwf_fp32 --mpi 2.2.2.2 --grid 48.48.48.48 --accelerator-threads 8 --shm 2048 $OPT > dwf.48.48.48.48.4node.opt0

OPT="--comms-overlap --comms-concurrent --shm-mpi 1"
srun ./benchmarks/Benchmark_dwf_fp32 --mpi 2.2.2.2 --grid 64.64.64.64 --accelerator-threads 8 --shm 2048 $OPT > dwf.64.64.64.64.4node.opt1
srun ./benchmarks/Benchmark_dwf_fp32 --mpi 2.2.2.2 --grid 48.48.48.48 --accelerator-threads 8 --shm 2048 $OPT > dwf.48.48.48.48.4node.opt1
