#!/bin/bash
#BSUB -P LGT104
#BSUB -W 0:20
#BSUB -nnodes 16
#BSUB -J DWF


export OMP_NUM_THREADS=6
export PAMI_IBV_ADAPTER_AFFINITY=1
export PAMI_ENABLE_STRIPING=1

DIR=.
source sourceme.sh

echo MPICH_SMP_SINGLE_COPY_MODE $MPICH_SMP_SINGLE_COPY_MODE

VOLS=( 32.32.32.16 32.32.32.64 64.32.32.64 64.32.64.64 64.64.64.64 64.64.64.128  64.64.64.256  64.64.64.512 128.64.64.64.512)
MPI=( 1.1.1.1      1.1.1.4     2.1.1.4         2.1.2.4     2.2.2.4      2.2.2.8      2.2.2.16      2.2.2.32 4.4.2.32 )
RANKS=(     1            4           8              16          32          64            128           256 1024)
NODES=(     1            1           2               4           8           16            32            64  128)
INTS=(      0            1           2               3           4            5             6             7    8)

for i in 5
do
    vol=${VOLS[$i]} 
    nodes=${NODES[$i]} 
    mpi=${MPI[$i]} 
    ranks=${RANKS[$i]} 

    JSRUN="jsrun --nrs $nodes -a4 -g4 -c42 -dpacked -b packed:10 --latency_priority gpu-cpu --smpiargs=-gpu"

    PARAMS=" --accelerator-threads 8 --grid $vol --mpi $mpi --comms-sequential --shm 2048 --shm-mpi 0"
    $JSRUN ./benchmarks/Benchmark_dwf_fp32 $PARAMS > run.v${vol}.n${nodes}.m${mpi}.seq.ker

    PARAMS=" --accelerator-threads 8 --grid $vol --mpi $mpi --comms-overlap --shm 2048 --shm-mpi 0"
    $JSRUN ./benchmarks/Benchmark_dwf_fp32 $PARAMS > run.v${vol}.n${nodes}.m${mpi}.over.ker

done

