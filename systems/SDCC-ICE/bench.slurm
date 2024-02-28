#!/bin/bash
#SBATCH --partition lqcd
#SBATCH --time=00:20:00
#SBATCH -A lqcdtest
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --qos lqcd

source sourceme.sh

export OMP_NUM_THREAD=24
#srun -N1 -n1 numactl -H > numa.txt
#srun -N1 -n1 lstopo ice-topo.pdf

cat << EOF > select_socket
#!/bin/bash
export NUM_MAP=(0 1)
export NUMA=\${NUMA_MAP[\$SLURM_LOCALID]}
exec \$*
EOF
chmod +x ./select_socket

#for vol in 8.8.8.16 8.8.8.32 8.8.8.64
#for vol in 8.8.16.16 8.8.16.32 8.8.16.64
for vol in 8.16.16.16 8.16.16.32 8.16.16.64 16.16.16.32 16.16.16.64 24.24.24.64 32.32.32.32
do
srun --cpu-bind=ldoms -N1 -n2 ./select_socket ./benchmarks/Benchmark_dwf_fp32 --mpi 1.1.1.2 --grid $vol --dslash-asm > $vol.2socket.out
srun --cpu-bind=ldoms -N1 -n1 ./select_socket ./benchmarks/Benchmark_dwf_fp32 --mpi 1.1.1.1 --grid $vol --dslash-asm > $vol.1socket.out
done

