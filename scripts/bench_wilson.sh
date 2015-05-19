for omp in 1 2 4
do
echo > wilson.t$omp
for vol in 4.4.4.4 4.4.4.8 4.4.8.8  4.8.8.8  8.8.8.8   8.8.8.16 8.8.16.16  8.16.16.16
do   
perf=` ./benchmarks/Grid_wilson --grid $vol --omp $omp  | grep mflop | awk '{print $3}'`
echo $vol $perf >> wilson.t$omp
done
done