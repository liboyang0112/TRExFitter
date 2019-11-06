#!/bin/bash
diff -w -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "png file FitExample" LOG_f test/logs/LOG_f && diff -w FitExample/Fits/FitExample.txt test/FitExample/Fits/FitExample.txt
