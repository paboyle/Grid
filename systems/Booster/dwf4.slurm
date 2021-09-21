#!/bin/sh
#SBATCH --account=gm2dwf
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=12
#SBATCH --time=2:00:00
#SBATCH --partition=develbooster
#SBATCH --gres=gpu:4

export OMP_NUM_THREADS=4
export OMPI_MCA_btl=^uct,openib
export UCX_TLS=gdr_copy,rc,rc_x,sm,cuda_copy,cuda_ipc
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=16384
export UCX_IB_GPU_DIRECT_RDMA=yes
export UCX_MEMTYPE_CACHE=n

OPT="--comms-overlap --comms-concurrent"

srun -N 4  -n $SLURM_NTASKS \
	./benchmarks/Benchmark_dwf_fp32 \
	$OPT \
	--mpi 2.2.2.2 \
	--accelerator-threads 8 \
	--grid 64.64.64.64 \
	--shm 2048 > dwf.4node.perf


srun -N 4  -n $SLURM_NTASKS \
	./benchmarks/Benchmark_comms_host_device \
	--mpi 2.2.2.2 \
	--accelerator-threads 8 \
	--grid 64.64.64.64 \
	--shm 2048 > comms.4node.perf




