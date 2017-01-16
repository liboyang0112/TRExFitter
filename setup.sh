# Setup ROOT and gcc
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
localSetupROOT 6.04.14-x86_64-slc6-gcc49-opt --quiet

if [ "${ROOTSYS}" == "" ]; then
   echo -e "\033[41;1;37m Error initializing ROOT. ROOT is not set up. Please check. \033[0m"
else
   echo -e "\033[42;1;37m ROOT has been set to: *${ROOTSYS}* \033[0m"
fi

# Michele
alias macro="root -l -b -q"
