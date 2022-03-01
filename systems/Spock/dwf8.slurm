#!/bin/bash
# Begin LSF Directives
#SBATCH -A LGT104
#SBATCH -t 3:00:00
#SBATCH -p ecp
#SBATCH -J DWF
#SBATCH -o DWF.%J
#SBATCH -e DWF.%J
#SBATCH -N 2
#SBATCH -n 8

DIR=.
module list
export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_SMP_SINGLE_COPY_MODE=CMA

export OMP_NUM_THREADS=8

AT=8
echo MPICH_SMP_SINGLE_COPY_MODE $MPICH_SMP_SINGLE_COPY_MODE
PARAMS=" --accelerator-threads ${AT} --grid 16.16.16.48 --mpi 1.2.2.2 --comms-overlap --shm 2048 --shm-mpi 0"
srun -N2 -n8 --label -c$OMP_NUM_THREADS --gpus-per-task=1 ./mpiwrapper.sh ./HMC/Mobius2p1f_DD_RHMC $PARAMS

