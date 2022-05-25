#!/bin/bash
# Begin LSF Directives
#SBATCH -A LGT104
#SBATCH -t 01:00:00
##SBATCH -U openmpThu
#SBATCH -J DWF
#SBATCH -o DWF.%J
#SBATCH -e DWF.%J
#SBATCH -N 8
#SBATCH -n 64
#SBATCH --exclusive
#SBATCH --gpu-bind=map_gpu:0,1,2,3,7,6,5,4

DIR=.
source setup.sh

export MPICH_OFI_NIC_POLICY=GPU
export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_SMP_SINGLE_COPY_MODE=XPMEM
#export MPICH_SMP_SINGLE_COPY_MODE=CMA
#export MPICH_SMP_SINGLE_COPY_MODE=NONE
export OMP_NUM_THREADS=1

echo MPICH_SMP_SINGLE_COPY_MODE $MPICH_SMP_SINGLE_COPY_MODE

for vol in 64.64.64.256 64.64.64.128 32.32.32.256  32.32.32.128
do
PARAMS=" --accelerator-threads 8 --grid $vol --mpi 2.2.2.8 --comms-overlap --shm 2048 --shm-mpi 1"
echo $PARAMS
srun --gpus-per-task 1 -N8 -n64  ./benchmarks/Benchmark_dwf_fp32 $PARAMS > dwf.${vol}.8node.shm-mpi1
done

PARAMS=" --accelerator-threads 8 --grid 64.64.64.32 --mpi 2.2.2.8 --comms-overlap --shm 2048 --shm-mpi 1"
echo $PARAMS
srun --gpus-per-task 1 -N8 -n64  ./benchmarks/Benchmark_ITT $PARAMS > itt.8node

PARAMS=" --accelerator-threads 8 --grid 64.64.64.32 --mpi 2.2.2.8 --comms-overlap --shm 2048 --shm-mpi 0"
echo $PARAMS
srun --gpus-per-task 1 -N8 -n64  ./benchmarks/Benchmark_ITT $PARAMS > itt.8node_shm0

