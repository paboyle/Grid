#!/bin/bash
##SBATCH -p PVC-SPR-QZEH 
##SBATCH -p PVC-ICX-QZNW

#SBATCH -p QZ1J-ICX-PVC

source /nfs/site/home/paboylex/ATS/GridNew/Grid/systems/PVC-nightly/setup.sh

export NT=16

# export IGC_EnableLSCFenceUGMBeforeEOT=0
# export SYCL_PROGRAM_COMPILE_OPTIONS="-ze-opt-large-register-file=False"
#export IGC_ShaderDumpEnable=1 
#export IGC_DumpToCurrentDir=1
export I_MPI_OFFLOAD=1
export I_MPI_OFFLOAD_TOPOLIB=level_zero
export I_MPI_OFFLOAD_DOMAIN_SIZE=-1
export SYCL_DEVICE_FILTER=gpu,level_zero
export I_MPI_OFFLOAD_CELL=tile
export EnableImplicitScaling=0
export EnableWalkerPartition=0
export SYCL_PI_LEVEL_ZERO_DEVICE_SCOPE_EVENTS=1
export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=1
export SYCL_PI_LEVEL_ZERO_USE_COPY_ENGINE=0

for i in 0 
do
mpiexec -launcher ssh -n 2 -host localhost  ./wrap4gpu.sh ./Benchmark_dwf_fp32 --mpi 1.1.1.2 --grid 32.32.32.64 --accelerator-threads $NT  --shm-mpi 1  --device-mem 32768
mpiexec -launcher ssh -n 2 -host localhost  ./wrap4gpu.sh ./Benchmark_dwf_fp32 --mpi 2.1.1.1 --grid 64.32.32.32 --accelerator-threads $NT  --shm-mpi 1  --device-mem 32768
done
#mpiexec -launcher ssh -n 2 -host localhost  ./wrap4gpu.sh ./Benchmark_halo --mpi 1.1.1.2 --grid 32.32.32.64 --accelerator-threads $NT  --shm-mpi 1 > halo.2tile.1x2.log
#mpiexec -launcher ssh -n 2 -host localhost  ./wrap4gpu.sh ./Benchmark_halo --mpi 2.1.1.1 --grid 64.32.32.32 --accelerator-threads $NT  --shm-mpi 1 > halo.2tile.2x1.log


