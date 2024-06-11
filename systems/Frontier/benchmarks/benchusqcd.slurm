#!/bin/bash -l
#SBATCH --job-name=bench
##SBATCH --partition=small-g
##SBATCH -q debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-node=8
#SBATCH --time=00:30:00
#SBATCH --account=phy157_dwf
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0

cat << EOF > select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3 7 6 5 4)
export NUMA_MAP=(3 3 1 1 2 2 0 0)
export GPU=\${GPU_MAP[\$SLURM_LOCALID]}
export NUMA=\${NUMA_MAP[\$SLURM_LOCALID]}
export HIP_VISIBLE_DEVICES=\$GPU
unset ROCR_VISIBLE_DEVICES
echo RANK \$SLURM_LOCALID using GPU \$GPU    
exec numactl -m \$NUMA -N \$NUMA \$*
EOF

chmod +x ./select_gpu

root=$HOME/Frontier/Grid/systems/Frontier/
source ${root}/sourceme.sh

export OMP_NUM_THREADS=7
export MPICH_GPU_SUPPORT_ENABLED=1
#export MPICH_SMP_SINGLE_COPY_MODE=XPMEM

srun ./select_gpu ./Benchmark_usqcd --grid 32.32.32.32 --mpi 1.2.2.2 --accelerator-threads 8 --comms-overlap --shm 4096 --shm-mpi 0 --grid $vol  > Benchmark_usqcd.log


