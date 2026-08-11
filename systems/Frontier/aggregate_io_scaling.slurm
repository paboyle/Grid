#!/bin/bash -l

# Node count scan.  Everything that could confound the comparison is held
# fixed and only the number of Lustre clients changes:
#
#   local volume   8.8.32.128        = 262144 sites x 576 B = 151 MB / rank
#   aggregation    k=2, row of 16, 8 extents of 18 MB        (identical at
#                                                             every node count)
#   old path       contiguous run = lLattice[0]*576 = 4.6 KB (identical too)
#
# So the ratio of contiguous extent between the two layout preserving paths is
# a constant 4000x across the whole scan, and any change in their relative
# bandwidth is a property of the client count alone.  That is the variable the
# 4 node runs could not move.
#
# Trim the tail of the run list if you want a smaller allocation; the
# geometries are independent.

#SBATCH --job-name=ioScale
#SBATCH --nodes=128
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-node=8
#SBATCH --time=02:00:00
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
unset ROCR_VISIBLE_DEVICES
export HIP_VISIBLE_DEVICES=\$GPU

rank=\$SLURM_PROCID
CMD="numactl -m \$NUMA -N \$NUMA  \$*"
exec \$CMD
EOF

chmod +x ./select_gpu
SELECT=$SLURM_SUBMIT_DIR/select_gpu

root=$HOME/ParallelIO/systems/Frontier/
source ${root}/sourceme-rocm7.2.sh

export OMP_NUM_THREADS=7
export MPICH_GPU_SUPPORT_ENABLED=1

WORK=/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/aggscale.$SLURM_JOB_ID
mkdir -p $WORK
cd       $WORK

# The largest rung writes three files of 154 GB, so budget ~700 GB for the
# whole scan and check the quota before submitting.
#
# Do NOT stripe by default.  Measured on Frontier at 128 nodes / 1024 ranks,
# 154.6 GB record, best of 3 (MB/s):
#
#                     default layout    lfs setstripe -c -1 -S 8M
#   raw                    97673                81728
#   MPI-IO lexicographic   31675                77534
#   aggregate             189535                86063
#
# The two layout preserving paths want opposite things.  Wide striping is what
# rescues MPI-IO at scale, 32 -> 78 GB/s, which is presumably what careful
# striping has been buying in production.  It costs the aggregate path a factor
# of 2.2, and costs raw as well.  With the default layout MPI-IO plateaus at
# ~30 GB/s from 64 nodes upward -- the 20-40 GB/s seen in production -- while
# the aggregate path keeps scaling.  Uncomment only to reproduce that.
#lfs setstripe -c -1 -S 8M $WORK
lfs getstripe -d $WORK      # record the inherited layout; needed to read the numbers

BIN=$SLURM_SUBMIT_DIR/tests/IO/Test_aggregate_io
ARGS="--aggregate-target 4194304 --io-no-correctness --io-reps 3"

run () {   # run <nodes> <grid> <mpi> <comment>
  local nodes=$1 nranks=$(( $1 * 8 ))
  echo
  echo "==================================================================="
  echo "=== SCAN nodes=$nodes ranks=$nranks  grid=$2  mpi=$3   $4"
  echo "==================================================================="
  srun -N$nodes -n$nranks --ntasks-per-node=8 $SELECT $BIN --grid $2 --mpi $3 \
       --accelerator-threads 8 --shm 2048 --device-mem 32000 $ARGS
  echo "=== exit $?"
}

#            nodes  global lattice      decomposition   record
run    4     32.32.64.128    4.4.2.1   "4.8 GB"
run    8     32.32.64.256    4.4.2.2   "9.7 GB"
run   16     32.32.64.512    4.4.2.4   "19.3 GB"
run   32     32.32.128.512   4.4.4.4   "38.6 GB"
run   64     32.32.128.1024  4.4.4.8   "77.3 GB"
run  128     32.32.256.1024  4.4.8.8   "154.6 GB"

echo
echo "=== done.  Files left in $WORK"
echo "=== For cold read numbers, submit a SECOND job on a fresh allocation"
echo "=== with WORK pointed at this directory and --io-read-only added to"
echo "=== ARGS.  Reading back what this job wrote only measures page cache."
ls -l $WORK
