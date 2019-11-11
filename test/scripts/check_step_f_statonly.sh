#!/bin/bash
diff -w -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "png file FitExample" LOG_STATONLY_f test/logs/LOG_STATONLY_f && diff -w FitExampleStatOnly/Fits/FitExampleStatOnly_statOnly.txt test/FitExampleStatOnly/Fits/FitExampleStatOnly_statOnly.txt
