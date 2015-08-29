#!/bin/bash

LOG=$1
SWEEPS=`grep dH $LOG | wc -l`
SWEEPS=`expr $SWEEPS - 80`
echo
echo $SWEEPS thermalised sweeps
echo
plaq=`grep Plaq $LOG | tail -n $SWEEPS | awk '{ S=S+$10} END { print S/NR} ' `
echo Plaquette: $plaq
echo
edH=`grep dH $LOG | tail -n $SWEEPS | awk '{ S=S+exp(-$10)} END { print S/NR} '`
echo "<e-dH>: $edH"

TRAJ=`grep Acc $LOG | wc -l`
ACC=`grep Acc $LOG | grep ACCE | wc -l`
PACC=`expr  100 \* ${ACC} / ${TRAJ} `
echo
echo "Acceptance $PACC %  $ACC / $TRAJ "

grep Plaq $LOG | awk '{ print $10 }' | uniq > plaq.dat
grep dH $LOG | awk '{ print $10 }' > dH.dat
echo set yrange [-1:1] > plot.gnu
echo set terminal 'pdf' >> plot.gnu
echo "set output 'plaq.${LOG}.pdf'" >> plot.gnu
echo "plot 'plaq.dat' w l, 'dH.dat' w l " >> plot.gnu
echo
gnuplot plot.gnu >& gnu.errs
open plaq.${LOG}.pdf


