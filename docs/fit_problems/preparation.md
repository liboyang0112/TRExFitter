# Hands on preparation

For the purpose of this tutorial, be sure to download and compile the latest version of TRExFitter (assuming you already ran `setupATLAS` and `lsetup git`)

```bash
git clone ssh://git@gitlab.cern.ch:7999/TRExStats/TRExFitter.git
cd TRExFitter

# checkout the release used for the tutorial
git checkout TtHFitter-00-04-05

# initialise submodules
git submodule init # You need to this only the first time
git submodule update

# source ROOT, cmake
source setup.sh

# compile the code
mkdir build
cd build
cmake ../
make -j4
cd ..

# copy the tutorial config files
cp /eos/user/t/tdado/TRExFitterTutorial/Configs/* config/

# you can run the code easily like this
trex-fitter h config/<ConfigFileName>
```
