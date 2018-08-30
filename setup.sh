#!bin/sh

if [ "${BASH_SOURCE[0]}" != "" ]; then
    # This should work in bash.
    _src=${BASH_SOURCE[0]}
elif [ "${ZSH_NAME}" != "" ]; then
    # And this in zsh.
    _src=${(%):-%x}
elif [ "${1}" != "" ]; then
    # If none of the above works, we take it from the command line.
    _src="${1/setup.sh/}/setup.sh"
else
    echo -e "\033[1;31mERROR:\033[0m Could not determine the base directory of TRExFitter, i.e. where \"setup.sh\" is located."
    echo -e "\033[1;31mERROR:\033[0m Can you give it to the source script as additional argument?"
    echo -e "\033[1;31mERROR:\033[0m For example: source ../setup.sh .."
    return 1
fi

location="$(cd -P "$(dirname "${_src}")" && pwd)"
unset _src


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
if [ ! "$(ls -A ${TREXFITTER_HOME}/CommonSystSmoothingTool)" ]; then
  echo -e "\033[1;31mERROR:\033[0m CommonSystSmoothingTool directory does not exist or is empty. "
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule init' (first time use) in the base directory of TRExFitter"
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule update' in the base directory of TRExFitter"
  return
fi

# Check if the CommonStatTools code exists
if [ ! "$(ls -A ${TREXFITTER_HOME}/CommonStatTools)" ]; then
  echo -e "\033[1;31mERROR:\033[0m CommonStatTools directory does not exist or is empty. "
  echo -e "\033[1;31mERROR:\033[0m You need to type 'git submodule init' (first time use) in the base directory of TRExFitter"
  echo -e "\033[1;31mERROR:\033[0m Followed with 'git submodule update' in the base directory of TRExFitter"
  return
fi

echo "Setting up cmake with: lsetup cmake"
lsetup cmake

export PATH=${PATH}${PATH:+:}${TREXFITTER_HOME}/build/bin

echo -e "\nIf everything went smoothly, continue by compiling via cmake (in a dedicated build directory):"
echo -e "mkdir -p build && cd build && cmake .. && cmake --build ."