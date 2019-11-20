#!/bin/bash
diff -w -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple" LOG_NTUPLE_p test/logs/LOG_NTUPLE_p && for file in `ls FitExampleNtuple/Tables/*postFit.txt`; do diff -w $file test/$file; done
