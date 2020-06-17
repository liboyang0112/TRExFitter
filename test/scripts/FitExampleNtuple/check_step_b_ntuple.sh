#!/bin/bash
diff -w -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple" LOG_NTUPLE_b test/logs/FitExampleNtuple/LOG_NTUPLE_b
