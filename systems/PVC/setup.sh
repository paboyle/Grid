export https_proxy=http://proxy-chain.intel.com:911
export LD_LIBRARY_PATH=/nfs/site/home/azusayax/install/lib:$LD_LIBRARY_PATH

module load intel-release
source /opt/intel/oneapi/PVC_setup.sh
#source /opt/intel/oneapi/ATS_setup.sh
module load intel/mpich/pvc45.3
export PATH=~/ATS/pti-gpu/tools/onetrace/:$PATH

#clsh embargo-ci-neo-022845
#source /opt/intel/vtune_amplifier/amplxe-vars.sh
