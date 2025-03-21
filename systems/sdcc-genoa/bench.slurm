#!/bin/bash
#SBATCH --partition lqcd
#SBATCH --time=00:50:00
#SBATCH -A lqcdtest
#SBATCH -q lqcd
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -w genoahost001,genoahost003,genoahost050,genoahost054
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --qos lqcd

source sourceme.sh

export PLACES=(1:16:4 1:32:2 0:64:1);
export THR=(16 32 64)

for t in  2 
do
   
export OMP_NUM_THREADS=${THR[$t]}
export OMP_PLACES=${PLACES[$t]}
export thr=${THR[$t]}

#for vol in 24.24.24.24 32.32.32.32 48.48.48.96
for vol in 48.48.48.96
do
srun -N1 -n1 ./benchmarks/Benchmark_dwf_fp32 --mpi 1.1.1.1 --grid $vol --dslash-asm --shm 8192 > $vol.1node.thr$thr
done
#srun -N1 -n1 ./benchmarks/Benchmark_usqcd --mpi 1.1.1.1 --grid $vol > usqcd.1node.thr$thr
done

