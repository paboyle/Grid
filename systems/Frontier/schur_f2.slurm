#!/bin/bash -l

#SBATCH --job-name=schurF2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-node=8
#SBATCH --time=00:30:00
#SBATCH --account=phy157_dwf
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH -q debug

#####################################################################
# F2 (schur_recursive_inverse_plan.txt 4C): the full DenseCoarseMatrix
# glue on GPU with a REAL coarsening -- DWF MdagM + 0.5 shift, random
# aggregation nbasis=4, per-dim blocks {4,4,3,3}, coarse 3.3.4.4/rank
# (N = 576 * nranks).  The test binary FORCES DENSE_SCHUR=2 audit mode
# and DENSE_PANEL_BYTES=65536 internally (setenv in main), so both
# invert paths run and are compared on the same imported A; multi-chunk
# panel gathers are exercised through the glue.
#
# First GPU execution in anger of: GRID_SIMT-guarded CoarsenOperator
# population kernel; fresh stencil->dense import (fp32 + fp64) with the
# transposed-contraction fix; fp64 rank-major import; distributed fp64
# recursion; single terminal fp32 rounding; boss rocSOLVER single-GCD
# path as the audit oracle.
#
# PRE-REGISTERED (laptop values, expect same magnitudes): operator
# health herm-dev ~1e-15; import certificate ~1e-8; fp64 import
# certificate ~5e-7 few-ulp; AUDIT rel ~3e-7 (pure fp32 rounding);
# VERIFY ~1.6e-7; T6 ALL PASS.
#####################################################################

cat << EOF > select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3 7 6 5 4)
export NUMA_MAP=(3 3 1 1 2 2 0 0)
export GPU=\${GPU_MAP[\$SLURM_LOCALID]}
export NUMA=\${NUMA_MAP[\$SLURM_LOCALID]}
unset ROCR_VISIBLE_DEVICES
export HIP_VISIBLE_DEVICES=\$GPU
CMD="numactl -m \$NUMA -N \$NUMA  \$*"
exec \$CMD
EOF
chmod +x ./select_gpu

# NUMA_MAP is only valid at 8 ranks/node (see aggregate_io.slurm);
# below 8 use the GCD-only wrapper and leave placement to Slurm.
cat << EOF > select_gcd
#!/bin/bash
export GPU_MAP=(0 1 2 3 7 6 5 4)
export GPU=\${GPU_MAP[\$SLURM_LOCALID]}
unset ROCR_VISIBLE_DEVICES
export HIP_VISIBLE_DEVICES=\$GPU
exec \$*
EOF
chmod +x ./select_gcd
SELECT=$SLURM_SUBMIT_DIR/select_gpu
SELECT_NONUMA=$SLURM_SUBMIT_DIR/select_gcd

root=$HOME/ParallelIO/systems/Frontier/
source ${root}/sourceme-rocm7.2.sh

export OMP_NUM_THREADS=7
export MPICH_GPU_SUPPORT_ENABLED=1

BIN=$SLURM_SUBMIT_DIR/tests/debug/Test_schur_dense_coarse

run () {   # run <ranks-per-node> <grid> <mpi> <comment>
  local rpn=$1 gr=$2 mp=$3 note=$4 wrap=$SELECT
  if [ $rpn -ne 8 ]; then wrap=$SELECT_NONUMA; fi
  echo
  echo "==================================================================="
  echo "=== F2  n=$rpn  grid=$gr  mpi=$mp   $note"
  echo "==================================================================="
  srun -N1 -n$rpn --ntasks-per-node=$rpn $wrap $BIN --grid $gr --mpi $mp \
       --accelerator-threads 8 --shm 2048 --device-mem 32000
  echo "=== exit $?"
}

# Uniform local fine volume 12.12.12.12 per rank (t = 12 * nranks).
# n=2 covers the padded-Extract branch; n=8 the full node, N=4608.
# The n=3 uneven case is laptop-proven; commented to keep the debug
# slot short -- uncomment for the full sweep.
run 2 12.12.12.24 1.1.1.2 "2 GCD, padded Extract on t, N=1152"
#run 3 12.12.12.36 1.1.1.3 "3 GCD, UNEVEN, N=1728"
run 8 12.12.12.96 1.1.1.8 "full node, 8 GCD, N=4608"

echo
echo "=== F2 done"
