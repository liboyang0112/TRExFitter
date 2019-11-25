#!/bin/bash
diff -w -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "png file FitExampleMorphing" LOG_MORPH_f test/logs/FitExampleMorphing/LOG_MORPH_f && diff -w FitExampleMorphing/Fits/FitExampleMorphing.txt test/reference/FitExampleMorphing/Fits/FitExampleMorphing.txt
