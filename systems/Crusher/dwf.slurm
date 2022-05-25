#!/bin/bash
# Begin LSF Directives
#SBATCH -A LGT104
#SBATCH -t 01:00:00
##SBATCH -U openmpThu
#SBATCH -J DWF
#SBATCH -o DWF.%J
#SBATCH -e DWF.%J
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --exclusive
#SBATCH --gpu-bind=map_gpu:0,1,2,3,7,6,5,4

DIR=.
source sourceme.sh
export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_SMP_SINGLE_COPY_MODE=XPMEM
export OMP_NUM_THREADS=1

echo MPICH_SMP_SINGLE_COPY_MODE $MPICH_SMP_SINGLE_COPY_MODE

echo working directory
pwd

PARAMS=" --accelerator-threads 8 --grid 32.32.32.32 --mpi 1.1.1.1 --comms-sequential --shm 2048 --shm-mpi 0"
srun --gpus-per-task 1 -n1  ./benchmarks/Benchmark_dwf_fp32 $PARAMS

PARAMS=" --accelerator-threads 8 --grid 64.64.64.32 --mpi 2.2.2.1 --comms-sequential --shm 2048 --shm-mpi 0"
srun --gpus-per-task 1 -n8  ./benchmarks/Benchmark_dwf_fp32 $PARAMS


