# Setting up TRExFitter

## Obtaining the package

Assuming your machine uses CentOS7 and has access to cvmfs (e.g. lxplus) you should first do
```bash
setupATLAS
lsetup git
```

To obtain the code simply clone the repository (you can also use different protocol)
```bash
git clone https://:@gitlab.cern.ch:8443/TRExStats/TRExFitter.git
```

If it is your first time clonig the code, you also need to download the submodules using
```bash
git submodule init
git submodule update
```

Each time the submodule changes, you need to repeat the `git submodule update` command.

To checkout a specific tag, simply do
```bash
git checkout <tag>
```

## Compiling the code

Go to the repository (`cd TRExFitter`) and then source the `setup.sh` script
```bash
source setup.sh
```

This will setup the correct ROOT and cmake versions.
Note that the script also supports old `slc6` version, which you can setup using `source setup.sh slc6`

To compile the code, first create a build repository, `mkdir build` and then go to the directory (`cd build`).
To compile the code do
```bash
cmake ../
make -j4
```

Or you can simply use an alias: `trex-make`.
The binary executable will appear in `build/bin/` folder.
