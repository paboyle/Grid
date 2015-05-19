plot 'wilson.t1' u 2 w l t "AVX1-OMP=1"
replot 'wilson.t2' u 2 w l t "AVX1-OMP=2"
replot 'wilson.t4' u 2 w l t "AVX1-OMP=4"
set terminal 'pdf'
set output 'wilson_clang.pdf'
replot
quit
