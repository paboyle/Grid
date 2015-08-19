#!/bin/bash

LOG=$1
SWEEPS=$2
grep Plaq $LOG | tail -n $SWEEPS | awk '{ S=S+$10} END { print S/NR} ' 
grep dH $LOG | tail -n $SWEEPS | awk '{ S=S+exp(-$10)} END { print S/NR} '
