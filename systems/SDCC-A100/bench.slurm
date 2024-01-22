#!/bin/bash
#SBATCH --partition csi
#SBATCH --time=00:10:00
#SBATCH -A csigeneral
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --qos csi
#SBATCH --gres=gpu:4

source sourceme.sh

cat << EOF > select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3)
export GPU=\${GPU_MAP[\$SLURM_LOCALID]}
export CUDA_VISIBLE_DEVICES=\$GPU
unset ROCR_VISIBLE_DEVICES
echo RANK \$SLURM_LOCALID using GPU \$GPU    
exec \$*
EOF
chmod +x ./select_gpu


export OMP_NUM_THREADS=4
export OMPI_MCA_btl=^uct,openib
export UCX_TLS=cuda,gdr_copy,rc,rc_x,sm,cuda_copy,cuda_ipc
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=16384
export UCX_IB_GPU_DIRECT_RDMA=no
export UCX_MEMTYPE_CACHE=n

export OMP_NUM_THREAD=8
#srun -N1 -n1 nvidia-smi
#srun -N1 -n1 numactl -H > numa.txt
srun -N1 -n1 lstopo A100-topo.pdf

# 4.35 TF/s
#srun -N1 -n1 ./benchmarks/Benchmark_dwf_fp32 --mpi 1.1.1.1 --grid 16.32.32.32 --shm 2048 --shm-mpi 0  --accelerator-threads 16

srun -N1 -n4 ./select_gpu ./benchmarks/Benchmark_dwf_fp32 --mpi 1.1.2.2 --grid 32.32.64.64 --shm 2048 --shm-mpi 0  --accelerator-threads 16

