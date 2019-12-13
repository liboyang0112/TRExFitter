#!/bin/bash
diff -w -I "Real time"  -I "RooRealVar::" -I "mkdir" -I "libSM.so" -I "libASImage" -I "png file FitExampleMorphMultifit" LOG_MORPH_MULTI_f test/logs/FitExampleMorphMultifit/LOG_MORPH_MULTI_f && diff -w FitExampleMorphMultifit/Fits/FitExampleMorphMultifit.txt test/reference/FitExampleMorphMultifit/Fits/FitExampleMorphMultifit.txt
