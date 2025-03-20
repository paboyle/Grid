#!/bin/bash
#SBATCH --partition lqcd
#SBATCH --time=00:50:00
#SBATCH -A lqcdtest
#SBATCH -q lqcd
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH -w genoahost001,genoahost003,genoahost050,genoahost054
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=64
#SBATCH --qos lqcd

source sourceme.sh

export PLACES=(1:16:4 1:32:2 0:64:1);
export THR=(16 32 64)

nodes=2
mpi=1.1.1.2

for t in 2 
do
   
export OMP_NUM_THREADS=${THR[$t]}
export OMP_PLACES=${PLACES[$t]}
export thr=${THR[$t]}

#srun -N$nodes -n$nodes ./benchmarks/Benchmark_usqcd --mpi $mpi --grid 32.32.32.32 > usqcd.n$nodes.thr$thr

for vol in 64.64.64.128
do
srun -N$nodes -n$nodes ./benchmarks/Benchmark_dwf_fp32 --mpi $mpi --grid $vol --dslash-asm --comms-overlap --shm 8192 > $vol.n$nodes.overlap.thr$thr
done

done

