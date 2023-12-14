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
#export BIND="--cpu-bind=verbose,map_ldom:3,3,1,1,2,2,0,0"

DIR=.
source sourceme.sh
export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=16384
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_SMP_SINGLE_COPY_MODE=CMA
export OMP_NUM_THREADS=4
export MPICH_OFI_NIC_POLICY=GPU

#PARAMS=" --accelerator-threads 8 --grid 64.64.32.16 --mpi 2.2.2.1 --comms-sequential --shm 2048 --shm-mpi 1"
#srun --gpus-per-task 1 -n8 $BIND ./wrap.sh  ./benchmarks/Benchmark_dwf_fp32 $PARAMS

PARAMS=" --accelerator-threads 8 --grid 16.16.32.32 --mpi 1.1.1.1 --comms-sequential --shm 2048 --shm-mpi 1"
srun --gpus-per-task 1 -n1 $BIND ./wrap.sh  ./benchmarks/Benchmark_dwf_fp32 $PARAMS

PARAMS=" --accelerator-threads 8 --grid 32.16.32.32 --mpi 1.1.1.2 --comms-sequential --shm 2048 --shm-mpi 1"
srun --gpus-per-task 1 -n2 $BIND ./wrap.sh  ./benchmarks/Benchmark_dwf_fp32 $PARAMS

PARAMS=" --accelerator-threads 8 --grid 32.32.32.64 --mpi 1.2.2.2 --comms-sequential --shm 2048 --shm-mpi 1"
srun --gpus-per-task 1 -n8 $BIND ./wrap.sh  ./benchmarks/Benchmark_dwf_fp32 $PARAMS


