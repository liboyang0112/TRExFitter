#!/bin/bash
diff -w -I "Real time"  -I "RooRealVar::" -I "mkdir" LOG_f test/logs/LOG_f && diff -w FitExample/Fits/FitExample.txt test/FitExample/Fits/FitExample.txt
