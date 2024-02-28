#!/bin/bash -l
#SBATCH --job-name=bench_lehner
#SBATCH --partition=small-g
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-node=8
#SBATCH --time=00:10:00
#SBATCH --account=project_465000546
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0

CPU_BIND="map_cpu:48,56,32,40,16,24,1,8"
echo $CPU_BIND

cat << EOF > select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3 4 5 6 7)
export GPU=\${GPU_MAP[\$SLURM_LOCALID]}
export HIP_VISIBLE_DEVICES=\$GPU
unset ROCR_VISIBLE_DEVICES
echo RANK \$SLURM_LOCALID using GPU \$GPU    
exec \$*
EOF

chmod +x ./select_gpu

root=/scratch/project_465000546/boylepet/Grid/systems/Lumi
source ${root}/sourceme.sh

export OMP_NUM_THREADS=7
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_SMP_SINGLE_COPY_MODE=XPMEM

for vol in 16.16.16.64 32.32.32.64  32.32.32.128
do
srun --cpu-bind=${CPU_BIND} ./select_gpu ./Benchmark_dwf_fp32 --mpi 2.2.2.2 --accelerator-threads 8 --comms-overlap --shm 2048 --shm-mpi 0 --grid $vol  > log.shm0.ov.$vol
#srun --cpu-bind=${CPU_BIND} ./select_gpu ./Benchmark_dwf_fp32 --mpi 2.2.2.2 --accelerator-threads 8 --comms-overlap --shm 2048 --shm-mpi 1 --grid $vol  > log.shm1.ov.$vol

srun --cpu-bind=${CPU_BIND} ./select_gpu ./Benchmark_dwf_fp32 --mpi 2.2.2.2 --accelerator-threads 8 --comms-sequential --shm 2048 --shm-mpi 0 --grid $vol  > log.shm0.seq.$vol
#srun --cpu-bind=${CPU_BIND} ./select_gpu ./Benchmark_dwf_fp32 --mpi 2.2.2.2 --accelerator-threads 8 --comms-sequential --shm 2048 --shm-mpi 1 --grid $vol > log.shm1.seq.$vol
done

