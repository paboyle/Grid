module load cuda/12.3 
module load ucx/1.15.0-cuda12.3  
module load openmpi/4.1.5-cuda12.3

export PREFIX=/home/tc002/tc002/shared/env/prefix/
export LD_LIBRARY_PATH=$PREFIX/lib/:$LD_LIBRARY_PATH
unset SBATCH_EXPORT

