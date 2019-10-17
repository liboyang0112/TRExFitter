#!/bin/bash

echo ""
echo "#########################################################################"
echo "TRExFitter validation logfiles"
echo "#########################################################################"
echo ""
echo "You are about to change the logfiles used online to validate the results."
echo "Note that ANY difference with the current logfiles must be justified in the"
echo "merge request form to make the bookmark of changes easier."
echo ""
echo "Are you sure you want to continue ? [y/n]"
read -n 1 ok_to_continue
echo ""
if [[ "${ok_to_continue}" != "y" ]]; then
  echo "Stopping the execution of the script !"
  return 0
fi
echo ""

echo "The generation of the new log files must be done from a *fresh terminal*, "
echo "removing the temporary compilation files ROOT is generating and *without "
echo "any setup done*. If you didn't do so, I will detect it and complain."
echo ""

##
## Check temporary cpmpilation files ... and deletes them if user is ok
##
for temp_file in `ls *_C_* *_c_*`; do
  echo "-> I found this file: $temp_file and will need to delete it. Ok ? [y/n]";
  read -n 1 ok_to_delete;
  echo "";
  if [[ "${ok_to_delete}" == "y" ]]; then
    rm $temp_file;
  fi
done

##
## Checks if ROOT is setup up by any other way as the setup.sh of the folder
##
which_root=`which root`
which_root=""
if [[ "$which_root" != "" ]]; then
  echo "Looks like ROOT has been setup :( Stopping the execution of the script now !"
  return 0;
fi

##
## Sources ROOT with THE good way
##
echo ""
source setup.sh #setting up ROOT as defined in the package
which_root=`which root`
which_root="toto"#FIXME
if [[ "$which_root" == "" ]]; then
  echo "For some reasons, the setup failed (access to /cvmfs/ ?) Stopping the "
  echo "execution of the script now !"
  return 0;
fi

##
## Produces the input ROOT file used to generate the fit
##
echo ""
echo "Now, I need to generate the input file by executing CreateHistograms.C"
root -l -b -q CreateHistograms.C && echo "Done !"
if [[ ! -d ExampleInputs/ ]]; then
  echo "!!ERROR!! The ExampleInputs folder cannot be found. I don't really know why ?"
  return 0
fi
echo ""

##
## Compiling the code
##
echo ""
echo "Compiling the code. Cleaning it first, and then recompile."
rm -rf build/
mkdir build && cd build/
cmake ../
cmake --build ./
cd ..
if [[ ! -f build/bin/trex-fitter ]]; then
  echo "!!ERROR!! The binary file is not found. Need to investigate !!"
  return 0
fi
echo ""

##
## Now, actually runs the thing
##
echo ""
echo "Things seem to be in order now. I am about to produce the logfiles. Note that "
echo "you still can abort the process at any moment."
echo ""
for step in h d w f l s d p ; do
  echo "==> $step step ongoing"
  ./build/bin/trex-fitter $step config/myFit.config >& LOG_$step
  cat LOG_$step | grep -v "TRExFitter" >& test/logs/LOG_$step
  rm -f LOG_$step
done

##
## Making a git status and asks if the files have to be added
##
echo ""
echo "Logfiles have been produced. So far, nothing has been added to git. This is "
echo "your responsability to do so. A few advices: "
echo "  - use git diff <path to the file> to see what has changed."
echo "  - if you understand ALL changes, do git add <path to the file>"
echo ""
echo "In case the post-fit yields/tables/systematics are expected to be different,"
echo "you will need to update also the corresponding files in test/FitExample."
echo ""
echo "When doing git commit: "
echo "  - please mention which test files have been updated and why in the commit"
echo "    message"
echo ""
echo "Push the change to your branch on the git server, and check the result of the"
echo "build, here: https://gitlab.cern.ch/TRExStats/TRExFitter/pipelines"
echo ""
