#!/bin/bash
diff -w -I "libSM.so" -I "libASImage" -I "png file FitExample" LOG_p test/logs/LOG_p && for file in `ls FitExample/Tables/*postFit.txt`; do diff -w $file test/$file; done
