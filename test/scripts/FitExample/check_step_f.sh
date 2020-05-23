#!/bin/bash
diff -w -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "NumericIntegration" -I "png file FitExample" LOG_f test/logs/FitExample/LOG_f && diff -w FitExample/Fits/FitExample.txt test/reference/FitExample/Fits/FitExample.txt
