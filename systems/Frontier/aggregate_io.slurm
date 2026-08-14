#!/bin/bash -l

#SBATCH --job-name=runIO
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-node=8
#SBATCH --time=02:00:00
#SBATCH --account=phy157_dwf
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH -q debug

cat << EOF > select_gpu
#!/bin/bash
export GPU_MAP=(0 1 2 3 7 6 5 4)
export NUMA_MAP=(3 3 1 1 2 2 0 0)
export GPU=\${GPU_MAP[\$SLURM_LOCALID]}
export NUMA=\${NUMA_MAP[\$SLURM_LOCALID]}
unset ROCR_VISIBLE_DEVICES
export HIP_VISIBLE_DEVICES=\$GPU

rank=\$SLURM_PROCID
CMD="numactl -m \$NUMA -N \$NUMA  \$*"
exec \$CMD
EOF

chmod +x ./select_gpu

# GPU_MAP/NUMA_MAP above describe a full node.  At fewer than 8 ranks per node
# Slurm hands each task cpus-per-task cores block-wise from core 0, so the NUMA
# domain NUMA_MAP names is not in the task's allowed cpu mask and numactl fails
# with EINVAL before Grid starts.  The GPU choice stays valid; only the memory
# binding does not.  For the small correctness runs the binding is irrelevant,
# so use a wrapper that picks the GCD and leaves placement to Slurm.
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

# Test_aggregate_io writes ref.bin / agg.bin / trunc.bin with *relative* paths,
# so they land in the cwd of each rank, which is the submission directory unless
# we move.  Home is NFS; the whole point of the exercise is Lustre.
WORK=/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/aggio.$SLURM_JOB_ID
mkdir -p $WORK
cd       $WORK

# Striping has to be set on the directory before the files are created, and
# all three paths inherit it, so the A/B comparison stays fair either way.
# Left unset: at 128 nodes the aggregate path is 2.2x faster on the default
# layout than on -c -1, while MPI-IO is 2.4x slower on it.  See the table in
# aggregate_io_scaling.slurm.  Set it only to reproduce that interaction.
#lfs setstripe -c 32 -S 8M $WORK
lfs getstripe -d $WORK      # record the inherited layout; needed to read the numbers

BIN=$SLURM_SUBMIT_DIR/tests/IO/Test_aggregate_io

# ROMIO / Cray MPICH collective buffering report for the reference path.
# Verbose; switch on when you want to see what the old path is doing.
# export MPICH_MPIIO_STATS=1
# export MPICH_MPIIO_TIMERS=1

# Always state ntasks-per-node explicitly.  Letting srun derive it from -n
# gives an uneven packing (3,2,2,2 for -n9 over 4 nodes) which no fixed
# localid -> GCD map can describe.  Ranks per node must not exceed 8.
run () {   # run <nodes> <ranks-per-node> <grid> <mpi> <comment> [extra args...]
  local nodes=$1 rpn=$2 gr=$3 mp=$4 note=$5 wrap=$SELECT
  local nranks=$(( $1 * $2 ))
  shift 5
  if [ $rpn -ne 8 ]; then wrap=$SELECT_NONUMA; fi   # NUMA_MAP only valid at 8
  echo
  echo "==================================================================="
  echo "=== N=$nodes x $rpn = $nranks  grid=$gr  mpi=$mp   $note"
  echo "=== extra: $@"
  echo "==================================================================="
  srun -N$nodes -n$nranks --ntasks-per-node=$rpn $wrap $BIN --grid $gr --mpi $mp \
       --accelerator-threads 8 --shm 2048 --device-mem 32000 "$@"
  echo "=== exit $?"
}

#####################################################################
# Phase 0.  Non-uniform AllToAllV, correctness only.
#
# The variable-count branch of CartesianCommunicator::AllToAllV only runs
# when the send counts differ between destinations, which needs an odd
# process factor in an un-split (fast) dimension.  It does NOT need a toy
# volume -- see phase 1b, which hits it at 226 MB/rank.  The one small run
# kept here is the only way to reach the non-uniform case at k=1 with many
# extents per rank; it is a logic corner, carries --io-reps 0, and reports
# no bandwidth.  2.2.4.4 per rank is 9 KB and means nothing for I/O.
#####################################################################
run 3 3 6.6.4.4  3.3.1.1  "non-uniform at all four targets, k=1 and k=3" --io-reps 0

#####################################################################
# Phase 1.  Medium volume, 4.8 GB record, 151 MB/rank.
# Correctness plus performance at the production 4MB target only; the
# byte-for-byte comparison reads the whole file on one rank, so the
# full four target sweep is not affordable above a few GB.
#####################################################################
run 4 8 32.32.64.128  4.4.2.1  "4.8 GB, correctness + perf" --aggregate-target 4194304

#####################################################################
# Phase 1b.  Non-uniform AllToAllV at a volume that means something.
# 24 ranks = 3 nodes x 8, so the NUMA map stays valid.
#   48.48.64.64 / 3.2.2.2  k=2 row=6  <=6 extents of 13.5 MB  NON-UNIFORM
# 5.4 GB record, 226 MB/rank.  Correctness and performance together, so
# the variable-count path is covered where the bandwidth is also real.
#####################################################################
run 3 8 48.48.64.64  3.2.2.2  "5.4 GB, NON-UNIFORM counts, correctness + perf" --aggregate-target 4194304

#####################################################################
# Phase 2.  Production volume, 19.3 GB record, 604 MB/rank.
# Correctness off.  Modelled extents at the 4MB target (fobj = 576 B):
#   1.1.4.8   k=2  row=1   16 extents of 36 MB   <- old path gives 36 KB runs
#   4.4.2.1   k=2  row=16   8 extents of 72 MB   <- old path gives  9 KB runs
#   2.2.2.4   k=2  row=4    8 extents
#   8.4.1.1   k=2  row=32   4 extents
#####################################################################
run 4 8 64.64.64.128  1.1.4.8  "19.3 GB, baseline decomposition" --aggregate-target 4194304 --io-no-correctness
run 4 8 64.64.64.128  4.4.2.1  "19.3 GB, row of 32"             --aggregate-target 4194304 --io-no-correctness
run 4 8 64.64.64.128  2.2.2.4  "19.3 GB, row of 4"              --aggregate-target 4194304 --io-no-correctness
run 4 8 64.64.64.128  8.4.1.1  "19.3 GB, row of 32, x heavy"    --aggregate-target 4194304 --io-no-correctness

# 10.9 GB with non-uniform counts, 453 MB/rank, 27 MB extents, 24 ranks.
run 3 8 48.96.64.64   3.4.2.1  "10.9 GB, NON-UNIFORM counts"    --aggregate-target 4194304 --io-no-correctness

#####################################################################
# Phase 3.  Large volume, 43.5 GB record, 1.36 GB/rank.
# Host side buffering is 3 copies of the local slice, so ~4 GB/rank
# and ~33 GB/node.  Two repetitions only.
#####################################################################
run 4 8 64.64.96.192  1.1.4.8  "43.5 GB, baseline"  --aggregate-target 4194304 --io-no-correctness --io-reps 2
run 4 8 64.64.96.192  4.4.2.1  "43.5 GB, row of 32" --aggregate-target 4194304 --io-no-correctness --io-reps 2

echo
echo "=== done.  Output left in $WORK"
ls -l $WORK
