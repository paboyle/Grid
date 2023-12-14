#!/bin/bash -l
#SBATCH --job-name=fthmc3ge
#SBATCH --partition=small-g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
##SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=8
#SBATCH --time=2:00:00
#SBATCH --account=project_465000546
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0


#sbatch --dependency=afterany:$SLURM_JOBID fthmc3gev.slurm

CPU_BIND="map_ldom:3,3,1,1,0,0,2,2"
MEM_BIND="map_mem:3,3,1,1,0,0,2,2"
echo $CPU_BIND

cat << EOF > ./select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3 4 5 6 7)
export NUMA_MAP=(3 3 1 1 0 0 2 2)
export GPU=\${GPU_MAP[\$SLURM_LOCALID]}
export NUM=\${NUMA_MAP[\$SLURM_LOCALID]}
#export HIP_VISIBLE_DEVICES=\$GPU
export ROCR_VISIBLE_DEVICES=\$GPU
echo RANK \$SLURM_LOCALID using GPU \$GPU    
echo NUMA \$SLURM_LOCALID using NUMA \${NUM}
echo numactl -m \$NUM -N \$NUM \$*
exec numactl -m \$NUM -N \$NUM \$*
EOF
cat ./select_gpu

chmod +x ./select_gpu

root=/scratch/project_465000546/boylepet/Grid/systems/Lumi
source ${root}/sourceme.sh

export OMP_NUM_THREADS=7
export MPICH_SMP_SINGLE_COPY_MODE=CMA
export MPICH_GPU_SUPPORT_ENABLED=1

#cfg=`ls -rt ckpoint_*lat* | tail -n 1  `
#traj="${cfg#*.}"
#cfg=`ls -rt ckpoint_*lat* | tail -n 1  `
traj=0

vol=32.32.32.64
mpi=1.2.2.2
PARAMS="--mpi $mpi --accelerator-threads 16 --comms-sequential --shm 2048 --shm-mpi 0 --grid $vol"
#HMCPARAMS="--StartingType CheckpointStart --StartingTrajectory $traj --Trajectories 200"
HMCPARAMS="--StartingType ColdStart --StartingTrajectory $traj --Trajectories 20"

srun ./select_gpu ../FTHMC2p1f_3GeV $HMCPARAMS $PARAMS

