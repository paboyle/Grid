#!/bin/bash
# Begin LSF Directives
#SBATCH -A LGT104
#SBATCH -t 01:00:00
##SBATCH -U openmpThu
#SBATCH -J DWF
#SBATCH -o DWF.%J
#SBATCH -e DWF.%J
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --gpu-bind=map_gpu:0,1

DIR=.
source setup.sh

export MPICH_OFI_NIC_POLICY=GPU
export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_NUM_THREADS=16

echo MPICH_SMP_SINGLE_COPY_MODE $MPICH_SMP_SINGLE_COPY_MODE

srun --gpus-per-task 1 -N1 -n2  ./benchmarks/Benchmark_dwf_fp32 --mpi 1.1.1.2 --grid 16.16.32.64 --shm-mpi 1 --shm 2048 --comms-sequential --accelerator-threads 8


