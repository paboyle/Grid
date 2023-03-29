#!/bin/bash
# Begin LSF Directives
#SBATCH -A LGT104
#SBATCH -t 01:00:00
##SBATCH -U openmpThu
#SBATCH -p ecp
#SBATCH -J comms
#SBATCH -o comms.%J
#SBATCH -e comms.%J
#SBATCH -N 1
#SBATCH -n 2

DIR=.
module list
export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
#export MPICH_SMP_SINGLE_COPY_MODE=XPMEM
#export MPICH_SMP_SINGLE_COPY_MODE=CMA
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export OMP_NUM_THREADS=8

AT=8
echo MPICH_SMP_SINGLE_COPY_MODE $MPICH_SMP_SINGLE_COPY_MODE
PARAMS=" --accelerator-threads ${AT} --grid 64.64.32.32 --mpi 2.1.1.1 "
srun -n2 --label -c$OMP_NUM_THREADS --gpus-per-task=1 ./mpiwrapper.sh ./benchmarks/Benchmark_comms_host_device $PARAMS

