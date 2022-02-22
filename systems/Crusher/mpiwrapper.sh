#!/bin/bash

lrank=$SLURM_LOCALID

export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID

echo "`hostname` - $lrank device=$ROCR_VISIBLE_DEVICES binding=$BINDING"

$*



