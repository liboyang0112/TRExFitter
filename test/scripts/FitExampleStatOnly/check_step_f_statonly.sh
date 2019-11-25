#!/bin/bash
diff -w -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "png file FitExampleStatOnly" LOG_STATONLY_f test/logs/FitExampleStatOnly/LOG_STATONLY_f && diff -w FitExampleStatOnly/Fits/FitExampleStatOnly_statOnly.txt test/reference/FitExampleStatOnly/Fits/FitExampleStatOnly_statOnly.txt
