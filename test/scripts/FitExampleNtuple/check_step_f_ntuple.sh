#!/bin/bash
diff -w -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple" LOG_NTUPLE_f test/logs/FitExampleNtuple/LOG_NTUPLE_f && diff -w FitExampleNtuple/Fits/FitExampleNtuple.txt test/reference/FitExampleNtuple/Fits/FitExampleNtuple.txt
