#!/bin/bash
diff -w -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple" LOG_NTUPLE_d test/logs/FitExampleNtuple/LOG_NTUPLE_d && for file in `ls FitExampleNtuple/Tables/*.txt`; do diff -w $file test/reference/FitExampleNtuple/$file; done
