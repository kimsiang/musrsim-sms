#! /bin/bash
echo ${1}
source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
cd /home/chencheng/musrsim-sms/build/run
pwd
/home/chencheng/musrsim-sms/build/musrSim 1000_${1}.mac ${1}