#!/bin/bash
#SBATCH -A m5294_g
#m3886_g
#SBATCH -C gpu
#SBATCH -q debug
#SBATCH -t 0:30:00
#SBATCH -c 32
#SBATCH -N 4
#SBATCH -n 16
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --exclusive
#SBATCH --gpu-bind=none

export SLURM_CPU_BIND="cores"
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_RDMA_ENABLED_CUDA=1
export MPICH_GPU_IPC_ENABLED=1
export MPICH_GPU_EAGER_REGISTER_HOST_MEM=0
export MPICH_GPU_NO_ASYNC_MEMCPY=0
#export MPICH_SMP_SINGLE_COPY_MODE=CMA

cat << EOF > select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3)
export NUMA_MAP=( 0 1 2 3 )
export GPU=\$SLURM_LOCALID
export NUMA=\$SLURM_LOCALID
export CUDA_VISIBLE_DEVICES=\$GPU
exec numactl -m \$NUMA -N \$NUMA \$*
EOF

chmod +x ./select_gpu

OPT="--comms-sequential --shm-mpi 0"
VOL=64.64.32.32
srun ./select_gpu ./benchmarks/Benchmark_usqcd --mpi 2.2.2.2 --grid $VOL --device-mem 16000 --accelerator-threads 8 --shm 2048 $OPT > usqcd.log


