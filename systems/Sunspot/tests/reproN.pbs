#!/bin/bash

#PBS -l select=32:system=sunspot,place=scatter
#PBS -A LatticeQCD_aesp_CNDA
#PBS -l walltime=02:00:00
#PBS -N reproN
#PBS -k doe

#export OMP_PROC_BIND=spread
#unset OMP_PLACES

module load oneapi/eng-compiler/2023.05.15.003
module load mpich/51.2/icc-all-deterministic-pmix-gpu

# 56 cores / 6 threads ~9
export OMP_NUM_THREADS=6
export MPIR_CVAR_CH4_OFI_ENABLE_GPU_PIPELINE=1
#export MPIR_CVAR_CH4_OFI_GPU_PIPELINE_D2H_ENGINE_TYPE=0
#export MPIR_CVAR_CH4_OFI_GPU_PIPELINE_H2D_ENGINE_TYPE=0
#export MPIR_CVAR_CH4_OFI_GPU_PIPELINE_BUFFER_SZ=1048576
#export MPIR_CVAR_CH4_OFI_GPU_PIPELINE_THRESHOLD=131072
#export MPIR_CVAR_CH4_OFI_GPU_PIPELINE_NUM_BUFFERS_PER_CHUNK=16
#export MPIR_CVAR_CH4_OFI_GPU_PIPELINE_MAX_NUM_BUFFERS=16
#export MPIR_CVAR_GPU_USE_IMMEDIATE_COMMAND_LIST=1

export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=1
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE=1
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE_FOR_D2D_COPY=1

export GRID_PRINT_ENTIRE_LOG=0
export GRID_CHECKSUM_RECV_BUF=1
export GRID_CHECKSUM_SEND_BUF=0

export MPICH_OFI_NIC_POLICY=GPU

export MPIR_CVAR_ALLREDUCE_DEVICE_COLLECTIVE=0
export MPIR_CVAR_REDUCE_DEVICE_COLLECTIVE=0
export MPIR_CVAR_ALLREDUCE_INTRA_ALGORITHM=recursive_doubling
unset MPIR_CVAR_CH4_COLL_SELECTION_TUNING_JSON_FILE
unset MPIR_CVAR_COLL_SELECTION_TUNING_JSON_FILE
unset MPIR_CVAR_CH4_POSIX_COLL_SELECTION_TUNING_JSON_FILE

cd $PBS_O_WORKDIR

NN=`cat $PBS_NODEFILE | wc -l`
echo $PBS_NODEFILE
cat $PBS_NODEFILE

echo $NN nodes in node file
for n in `eval echo {1..$NN}`
do

cd $PBS_O_WORKDIR

THIS_NODE=`head -n$n $PBS_NODEFILE | tail -n1 `
echo Node $n is $THIS_NODE

DIR=reproN.$PBS_JOBID/node-$n-$THIS_NODE

mkdir -p $DIR
cd $DIR

echo $THIS_NODE > nodefile

#CMD="mpiexec -np 12 -ppn 12  -envall --hostfile nodefile \
#	     ../../gpu_tile_compact.sh \
#	     ../../Test_dwf_mixedcg_prec --mpi 1.2.2.3 --grid 32.64.64.96 \
#		--shm-mpi 0 --shm 4096 --device-mem 32000 --accelerator-threads 32 --seconds 6000 --debug-stdout --log Message --comms-overlap"

CMD="mpiexec -np 12 -ppn 12  -envall --hostfile nodefile \
	     ../../gpu_tile_compact.sh \
	     ../../Test_dwf_mixedcg_prec --mpi 1.2.2.3 --grid 32.64.64.96 \
		--shm-mpi 1 --shm 4096 --device-mem 32000 --accelerator-threads 32 --seconds 6000 --debug-stdout --log Message --comms-overlap"

echo $CMD > command-line
env > environment
$CMD &

done

# Suspicious wait is allowing jobs to collide and knock out
#wait

sleep 6500

for n in ` eval echo {1..$NN} `
do

THIS_NODE=`head -n$n $PBS_NODEFILE | tail -n1 `
DIR=reproN.$PBS_JOBID/node-$n-$THIS_NODE

cd $DIR

grep Oops Grid.stderr.* > failures.$PBS_JOBID
rm core.*

done
