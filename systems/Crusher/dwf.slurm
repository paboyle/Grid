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
module list
export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
#export MPICH_SMP_SINGLE_COPY_MODE=XPMEM
export MPICH_SMP_SINGLE_COPY_MODE=NONE
#export MPICH_SMP_SINGLE_COPY_MODE=CMA
export OMP_NUM_THREADS=1

echo MPICH_SMP_SINGLE_COPY_MODE $MPICH_SMP_SINGLE_COPY_MODE

PARAMS=" --accelerator-threads 16 --grid 32.32.32.256 --mpi 1.1.1.8 --comms-overlap --shm 2048 --shm-mpi 0"
echo $PARAMS
srun --gpus-per-task 1 -n8  ./benchmarks/Benchmark_dwf_fp32 $PARAMS



