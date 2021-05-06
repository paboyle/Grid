#!/bin/bash

LOG=$1
SWEEPS=`grep dH.= $LOG | wc -l`
SWEEPS=`expr $SWEEPS - 100`
echo
echo $SWEEPS thermalised sweeps
echo
plaq=`grep Plaq $LOG | tail -n $SWEEPS | awk '{ S=S+$12} END { print S/NR} ' `
plaqe=`grep Plaq $LOG | tail -n $SWEEPS | awk '{ S=S+$12 ; SS=SS+$12*$12 } END { print sqrt( (SS/NR - S*S/NR/NR)/NR) } ' `
echo "Plaquette: $plaq (${plaqe})"
echo

grep  Plaq $LOG | tail -n $SWEEPS | awk '{ S=S+$12/20; if(NR%20==0){ print NR/20, " ", S; S=0;} } '  > plaq.binned

plaq=`cat plaq.binned  | awk '{ S=S+$2} END { print S/NR} ' `
plaqe=`cat plaq.binned | awk '{ S=S+$2 ; SS=SS+$2*$2 } END { print sqrt( (SS/NR - S*S/NR/NR)/NR) } ' `
echo "Binned Plaquette: $plaq (${plaqe})"
echo

dHv=`grep dH.= $LOG | tail -n $SWEEPS | awk '{ S=S+$16 ; SS=SS+$16*$16 } END { print sqrt(SS/NR) } ' `
edH=`grep dH.= $LOG | tail -n $SWEEPS | awk '{ S=S+exp(-$16)} END { print S/NR} '`
dedH=`grep dH.= $LOG | tail -n $SWEEPS | awk '{ S=S+exp(-$16); SS=SS+exp(-$16)*exp(-$16)} END { print sqrt( (SS/NR - S*S/NR/NR)/NR) } '`
echo "<e-dH>: $edH (${dedH})"
echo "<rms dH>: $dHv"

TRAJ=`grep Acc $LOG | wc -l`
ACC=`grep Acc $LOG | grep ACCE | wc -l`
PACC=`expr  100 \* ${ACC} / ${TRAJ} `
echo
echo "Acceptance $PACC %  $ACC / $TRAJ "

grep Plaq $LOG | awk '{ print $12 }' | uniq > plaq.dat
grep dH.= $LOG | awk '{ print $16 }' > dH.dat
echo set yrange [0.58:0.60] > plot.gnu
echo set terminal 'pdf' >> plot.gnu
echo "f(x) =0.588" >> plot.gnu
echo "set output 'plaq.${LOG}.pdf'" >> plot.gnu
echo "plot 'plaq.dat' w l, f(x) " >> plot.gnu
echo
gnuplot plot.gnu >& gnu.errs
open plaq.${LOG}.pdf


