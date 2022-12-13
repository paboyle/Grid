#!/bin/sh
##SBATCH -p PVC-SPR-QZEH 
##SBATCH -p PVC-ICX-QZNW
#SBATCH -p QZ1J-ICX-PVC
##SBATCH -p QZ1J-SPR-PVC-2C

source /nfs/site/home/paboylex/ATS/GridNew/Grid/systems/PVC-nightly/setup.sh

export NT=8

export I_MPI_OFFLOAD=1
export I_MPI_OFFLOAD_TOPOLIB=level_zero
export I_MPI_OFFLOAD_DOMAIN_SIZE=-1

# export IGC_EnableLSCFenceUGMBeforeEOT=0
# export SYCL_PROGRAM_COMPILE_OPTIONS="-ze-opt-large-register-file=False"
export SYCL_DEVICE_FILTER=gpu,level_zero
#export IGC_ShaderDumpEnable=1 
#export IGC_DumpToCurrentDir=1
export I_MPI_OFFLOAD_CELL=tile
export EnableImplicitScaling=0
export EnableWalkerPartition=0
export ZE_AFFINITY_MASK=0.0
mpiexec -launcher ssh -n 1 -host localhost  ./Benchmark_dwf_fp32 --mpi 1.1.1.1 --grid 32.32.32.32 --accelerator-threads $NT --comms-sequential --shm-mpi 1 --device-mem 32768

export ZE_AFFINITY_MASK=0
export I_MPI_OFFLOAD_CELL=device
export EnableImplicitScaling=1
export EnableWalkerPartition=1




















#mpiexec -launcher ssh -n 2 -host localhost  vtune -collect gpu-hotspots -knob gpu-sampling-interval=1 -data-limit=0 -r ./vtune_run4 -- ./wrap.sh ./Benchmark_dwf_fp32 --mpi 2.1.1.1 --grid 64.32.32.32 --accelerator-threads $NT --comms-overlap --shm-mpi 1

#mpiexec  -launcher ssh -n 1 -host localhost ./wrap.sh ./Benchmark_dwf_fp32 --mpi 1.1.1.1 --grid 64.32.32.32 --accelerator-threads $NT --comms-overlap --shm-mpi 1

#mpiexec  -launcher ssh -n 2 -host localhost ./wrap.sh ./Benchmark_dwf_fp32 --mpi 2.1.1.1 --grid 64.32.32.32 --accelerator-threads $NT --comms-sequential --shm-mpi 1

#mpiexec  -launcher ssh -n 2 -host localhost ./wrap.sh ./Benchmark_dwf_fp32 --mpi 2.1.1.1 --grid 64.32.32.32 --accelerator-threads $NT --comms-overlap --shm-mpi 1

#mpiexec  -launcher ssh -n 2 -host localhost ./wrap.sh ./Benchmark_dwf_fp32 --mpi 2.1.1.1 --grid 64.32.32.32 --accelerator-threads $NT --comms-sequential --shm-mpi 0

#mpirun -np 2 ./wrap.sh ./Benchmark_dwf_fp32 --mpi 1.1.1.2 --grid 16.32.32.64 --accelerator-threads $NT --comms-sequential --shm-mpi 0
#mpirun -np 2 ./wrap.sh ./Benchmark_dwf_fp32 --mpi 1.1.1.2 --grid 32.32.32.64 --accelerator-threads $NT --comms-sequential --shm-mpi 1

