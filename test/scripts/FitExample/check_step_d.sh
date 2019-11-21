#!/bin/bash
diff -w -I "libSM.so" -I "libASImage" -I "png file FitExample" LOG_d test/logs/FitExample/LOG_d && for file in `ls FitExample/Tables/*.txt`; do diff -w $file test/reference/$file; done
