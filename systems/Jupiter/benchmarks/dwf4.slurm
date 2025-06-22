#!/bin/sh
#SBATCH --account=jureap14
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=64
#SBATCH --time=2:00:00
#SBATCH --partition=booster
#SBATCH --gres=gpu:4

export OMP_NUM_THREADS=4
export OMPI_MCA_btl=^uct,openib
export UCX_TLS=gdr_copy,rc,rc_x,sm,cuda_copy,cuda_ipc
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=16384
export UCX_IB_GPU_DIRECT_RDMA=yes
export UCX_MEMTYPE_CACHE=n

OPT="--comms-overlap"

source ../sourceme.sh

cat << EOF > bind_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3)
export NUMA_MAP=(0 1 2 3)
export NIC_MAP=(0 1 2 3)
export GPU=\$SLURM_LOCALID
export NUMA=\$SLURM_LOCALID
export NIC=\$SLURM_LOCALID
export CUDA_VISIBLE_DEVICES=\$GPU
export UCX_NET_DEVICES=mlx5_\${NIC}:1

echo RANK \$SLURM_LOCALID using NUMA \$NUMA  GPU \$GPU NIC \$UCX_NET_DEVICES
exec numactl -m \$NUMA -N \$NUMA \$*
EOF

chmod +x ./bind_gpu

srun --cpu-bind=no -N 4 -n $SLURM_NTASKS \
        ./bind_gpu ./Benchmark_dwf_fp32 \
	$OPT \
	--mpi 2.2.2.2 \
	--accelerator-threads 8 \
	--grid 64.64.64.64 \
	--shm 2048 > dwf.4node.perf

srun --cpu-bind=no -N 4  -n $SLURM_NTASKS \
	./bind_gpu ./Benchmark_comms_host_device \
	--mpi 2.2.2.2 \
	--accelerator-threads 8 \
	--grid 32.32.64.64 \
	--shm 2048 > comms.4node.perf




