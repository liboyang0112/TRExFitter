#!bin/bash
location=$1

# Setup ROOT and gcc
# added back by Michele
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
# removed
# source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
# localSetupROOT 6.04.14-x86_64-slc6-gcc49-opt --quiet
localSetupROOT 6.10.06-x86_64-slc6-gcc62-opt --quiet

if [ "${ROOTSYS}" = "" ]; then
   echo -e "\033[41;1;37m Error initializing ROOT. ROOT is not set up. Please check. \033[0m"
else
   echo -e "\033[42;1;37m ROOT has been set to: *${ROOTSYS}* \033[0m"
fi

alias macro="root -l -b -q"

if [[ "$location" != "" ]]
then
  export PATH=$PATH:$location
  # to be able to point to the confg schema
  export TREXFITTER_HOME=$location
else
  export PATH=$PATH:`pwd`
  # to be able to point to the confg schema
  export TREXFITTER_HOME=`pwd`
fi

if [ ! -f $TREXFITTER_HOME/logo.txt ]; then
  echo -e "\033[1;31mWARNING:\033[0m \$TREXFITTER_HOME environmental variable not set properly"
  echo -e "\033[1;31mWARNING:\033[0m call this script with the path to the TRExFitter directory as an additional argument"
fi

# Check if the CommomSmoothing code exists
if [ ! "$(ls -A CommonSystSmoothingTool)" ]; then
  echo -e "\033[1;31mERROR:\033[0m CommonSystSmoothingTool directory does not exist or is empty. "
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule update' "
  return
fi

# Check if the CommonStatTools code exists
if [ ! "$(ls -A CommonStatTools)" ]; then
  echo -e "\033[1;31mERROR:\033[0m CommonStatTools directory does not exist or is empty. "
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule update' "
  return
fi

echo "Setting up cmake with: lsetup cmake"
lsetup cmake

