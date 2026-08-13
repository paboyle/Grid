#!/bin/bash -l

# Standalone MPI-only I/O reproducer on Frontier.  No Grid, no accelerator,
# so no GCD/NUMA wrapper is needed -- the point of the exercise is that this
# depends on nothing but an MPI installation and a filesystem.
#
#   mpicxx -O2 -std=c++11 io_mpi.cc -o io_mpi
#
# Weak scan: the local volume, and therefore the file view structure, is held
# identical at every rung and only the number of Lustre clients changes:
#
#   local volume  8.8.32.128 = 262144 sites x 576 B = 151 MB/rank
#   file view     32768 contiguous runs of 4608 B per rank, at every rung
#   aggregate     k=2, row of 16, 8 extents of 18 MB, at every rung
#
# so any change in the relative bandwidth of the two lexicographic paths is a
# property of the client count alone.
#
# The PERF lines are MiB/s (bytes/1024/1024/s), which is what BinaryIO.h
# computes for lastPerf.mbytesPerSecond and prints as "MB/s", so the two
# tools can be compared directly.  Grid's timed region is used here too:
# barrier, start, [plan build + exchange + I/O], barrier, stop, quoting the
# boss rank's stopwatch.  --reuse-plan hoists the plan build out, which is
# how to show it is not where the time goes; do not use it when comparing
# against Grid's own numbers.
#
# io_aurora.pbs runs 12 ranks per node, one per tile, because that is what
# that machine is.  The per rank local volume and the file view are the same
# there as here, but the record size and client count at a given NODE count
# are 1.5x.  See the header of that script before tabulating the two
# together.

#SBATCH --job-name=ioMPI
#SBATCH --nodes=128
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --time=02:00:00
#SBATCH --account=phy157_dwf
#SBATCH --exclusive
#SBATCH --mem=0

module load cce/21.0.0
module load cpe/26.03

WORK=/lustre/orion/phy157/proj-shared/phy157_dwf/paboyle/iompi.$SLURM_JOB_ID
mkdir -p $WORK
cd       $WORK

# Do NOT stripe by default.  Wide striping is what rescues the collective at
# scale and costs every other path a factor of 1.2-2.2; the default layout is
# what a user gets without knowing to ask.  Uncomment to reproduce that
# interaction, and record which one you ran.
#lfs setstripe -c -1 -S 8M $WORK
lfs getstripe -d $WORK

BIN=$SLURM_SUBMIT_DIR/io_mpi
ARGS="--target 4194304 --reps 3"

# ROMIO's own view of what it did.  Verbose, but the first thing anyone
# reading the report will ask for.
# export MPICH_MPIIO_STATS=1
# export MPICH_MPIIO_TIMERS=1

run () {   # run <nodes> <grid> <mpi> <comment> [extra args...]
  local nodes=$1 gr=$2 mp=$3 note=$4
  local nranks=$(( nodes * 8 ))
  shift 4
  echo
  echo "==================================================================="
  echo "=== nodes=$nodes ranks=$nranks grid=$gr mpi=$mp   $note"
  echo "=== extra: $@"
  echo "==================================================================="
  srun -N$nodes -n$nranks --ntasks-per-node=8 $BIN --grid $gr --mpi $mp $ARGS "$@"
  echo "=== exit $?"
}

#####################################################################
# Phase 0.  Correctness, including the non-uniform Alltoallv branch
# (odd process factor in an un-split dimension).  Small, and the
# whole-file crc32 is serial, so keep the volume down here.
#####################################################################
run 1 12.12.8.8   2.2.2.1  "correctness, uniform counts"     --reps 0 --serial-crc
run 3 24.12.8.8   3.2.2.2  "correctness, NON-uniform counts" --reps 0 --serial-crc
run 4 16.16.16.32 2.2.2.4  "correctness, non-zero offset"    --reps 0 --serial-crc --offset 1024

#####################################################################
# Phase 1.  Weak scan, 151 MB/rank.  Timing only.
#####################################################################
run    4     32.32.64.128    4.4.2.1   "4.8 GB"    --no-validate
run    8     32.32.64.256    4.4.2.2   "9.7 GB"    --no-validate
run   16     32.32.64.512    4.4.2.4   "19.3 GB"   --no-validate
run   32     32.32.128.512   4.4.4.4   "38.6 GB"   --no-validate
run   64     32.32.128.1024  4.4.4.8   "77.3 GB"   --no-validate
run  128     32.32.256.1024  4.4.8.8   "154.6 GB"  --no-validate

#####################################################################
# Phase 2.  Answer the two questions a reviewer will ask immediately.
#####################################################################
# Is the gap an artefact of measuring cache rather than the filesystem?
run  128     32.32.256.1024  4.4.8.8   "154.6 GB, durable"  --no-validate --fsync --drop-cache
# Does the collective recover if it is given the hints it wants?
run  128     32.32.256.1024  4.4.8.8   "154.6 GB, cb hints" --no-validate \
     --hints romio_cb_write=enable,romio_cb_read=enable,cb_nodes=128,cb_buffer_size=16777216
# Does the degenerate memory subarray matter?
run  128     32.32.256.1024  4.4.8.8   "154.6 GB, mem subarray" --no-validate --mem-subarray

echo
echo "=== done.  Files left in $WORK"
ls -l $WORK
