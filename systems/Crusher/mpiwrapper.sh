#!/bin/bash

lrank=$SLURM_LOCALID
lgpu=(0 1 2 3 7 6 5 4)

export ROCR_VISIBLE_DEVICES=${lgpu[$lrank]}

echo "`hostname` - $lrank device=$ROCR_VISIBLE_DEVICES "

$*



