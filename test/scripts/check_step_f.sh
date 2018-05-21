#!/bin/bash
diff -I "-- Real time"  -I "RooRealVar::" LOG_f test/logs/LOG_f && diff FitExample/Fits/FitExample.txt test/FitExample/Fits/FitExample.txt
