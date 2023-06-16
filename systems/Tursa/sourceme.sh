module load cuda/11.4.1  openmpi/4.1.1-cuda11.4.1  ucx/1.12.0-cuda11.4.1  
#module load cuda/11.4.1 openmpi/4.1.1 ucx/1.10.1
export PREFIX=/home/tc002/tc002/shared/env/prefix/
export LD_LIBRARY_PATH=$PREFIX/lib/:$LD_LIBRARY_PATH
unset SBATCH_EXPORT

