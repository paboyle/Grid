#!/bin/bash
#SBATCH -J dslash
#SBATCH -A tc002
#SBATCH -t 2:20:00
#SBATCH --exclusive
#SBATCH --nodes=16
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:4
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

export OMP_NUM_THREADS=4
export OMPI_MCA_btl=^uct,openib
export UCX_TLS=gdr_copy,rc,rc_x,sm,cuda_copy,cuda_ipc
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=16384
export UCX_IB_GPU_DIRECT_RDMA=yes
export UCX_MEMTYPE_CACHE=n
OPT="--comms-overlap --comms-concurrent"


mpirun -np $SLURM_NTASKS -x LD_LIBRARY_PATH  --bind-to none ./mpiwrapper.sh \
       	./benchmarks/Benchmark_dwf_fp32 \
	$OPT \
	--mpi 2.2.2.8 \
       	--accelerator-threads 8 \
	--grid 64.64.64.256 \
	--shm 2048 > dwf.16node.perf

