#!/bin/bash
diff -w -I "/afs/cern.ch/user/l/lvalery/" -I "/builds/" -I "Real time" -I "mkdir" -I "RooRealVar" -I "ibSM.so" -I "libASImage" -I "png file FitExample" LOG_r test/logs/FitExample/LOG_r && diff -w FitExample/Fits/NPRanking_SigXsecOverSM.txt test/reference/FitExample/Fits/NPRanking_SigXsecOverSM.txt
