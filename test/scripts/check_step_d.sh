#!/bin/bash
diff -w LOG_d test/logs/LOG_d && for file in `ls FitExample/Tables/*.txt`; do diff -w $file test/$file; done
