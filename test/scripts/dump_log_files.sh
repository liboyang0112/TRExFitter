#!/bin/bash
for step in d w f l s d p ; do
  ./myFit.exe $step config/myFit.config >& LOG_$step
  cat LOG_$step | grep -v "TRExFitter" >& test/logs/LOG_$step
done
