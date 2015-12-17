#!/bin/bash

LOG=$1
SWEEPS=`grep dH $LOG | wc -l`
SWEEPS=`expr $SWEEPS - 80`
echo
echo $SWEEPS thermalised sweeps
echo
plaq=`grep Plaq $LOG | tail -n $SWEEPS | awk '{ S=S+$10} END { print S/NR} ' `
plaqe=`grep Plaq $LOG | tail -n $SWEEPS | awk '{ S=S+$10 ; SS=SS+$10*$10 } END { print sqrt( (SS/NR - S*S/NR/NR)/NR) } ' `
echo "Plaquette: $plaq (${plaqe})"
echo

dHv=`grep dH $LOG | tail -n $SWEEPS | awk '{ S=S+$10 ; SS=SS+$10*$10 } END { print sqrt(SS/NR) } ' `
edH=`grep dH $LOG | tail -n $SWEEPS | awk '{ S=S+exp(-$10)} END { print S/NR} '`
echo "<e-dH>: $edH"
echo "<rms dH>: $dHv"

TRAJ=`grep Acc $LOG | wc -l`
ACC=`grep Acc $LOG | grep ACCE | wc -l`
PACC=`expr  100 \* ${ACC} / ${TRAJ} `
echo
echo "Acceptance $PACC %  $ACC / $TRAJ "

grep Plaq $LOG | awk '{ print $10 }' | uniq > plaq.dat
grep dH $LOG | awk '{ print $10 }' > dH.dat
echo set yrange [-0.2:1.0] > plot.gnu
echo set terminal 'pdf' >> plot.gnu
echo "set output 'plaq.${LOG}.pdf'" >> plot.gnu
echo "plot 'plaq.dat' w l, 'dH.dat' w l " >> plot.gnu
echo
gnuplot plot.gnu >& gnu.errs
open plaq.${LOG}.pdf


