#!/bin/bash -l

#SBATCH --job-name=schurF1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-node=8
#SBATCH --time=00:15:00
#SBATCH --account=phy157_dwf
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH -q debug

#####################################################################
# F1 (schur_recursive_inverse_plan.txt 4C): RecursiveSchurInverse on
# GPU + GPU-aware comms -- no lattice physics, communicator only.
# Laptop-proven semantics (Eigen oracle, ranks 1-4); this run
# adjudicates ONLY the device paths: rocSOLVER zgetri leaves, ComplexD
# strided gemmBatched merges, device panels, thrust-complex NegateCopy
# kernel, multi-GCD GlobalSumVector panels.
#
# PRE-REGISTERED: every stage PASS; T4 certificates ~1e-14..1e-10
# (eps64 growth model, as laptop); at P=1-equivalent leaves rocSOLVER
# != Eigen so |X-Xref| = exact 0 is NOT expected, <=1e-12 is healthy.
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

BIN=$SLURM_SUBMIT_DIR/tests/debug/Test_schur_inverse

run () {   # run <ranks-per-node> <grid> <mpi> <comment>
  local rpn=$1 gr=$2 mp=$3 note=$4 wrap=$SELECT
  if [ $rpn -ne 8 ]; then wrap=$SELECT_NONUMA; fi
  echo
  echo "==================================================================="
  echo "=== F1  n=$rpn  grid=$gr  mpi=$mp   $note"
  echo "==================================================================="
  srun -N1 -n$rpn --ntasks-per-node=$rpn $wrap $BIN --grid $gr --mpi $mp \
       --accelerator-threads 8 --shm 2048 --device-mem 32000
  echo "=== exit $?"
}

# Rank-count sweep: even splits, the uneven n=3 case, and the full node.
run 2 16.16.16.16 1.1.1.2 "even split, 2 GCD"
run 3 16.16.16.24 1.1.1.3 "UNEVEN split, 3 GCD"
run 4 16.16.16.16 1.1.1.4 "even split, 4 GCD"
run 8 16.16.16.16 1.1.1.8 "full node, 8 GCD, depth-3 tree"

echo
echo "=== F1 done"
