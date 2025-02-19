#!/bin/sh

export HIP_VISIBLE_DEVICES=$ROCR_VISIBLE_DEVICES
unset ROCR_VISIBLE_DEVICES

#rank=$SLURM_PROCID
#rocprof -d rocprof.$rank -o rocprof.$rank/results.rank$SLURM_PROCID.csv --sys-trace $@

$@
