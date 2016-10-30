#!/bin/bash
diff LOG_p test/logs/LOG_p && for file in `ls FitExample/Tables/*postFit.txt`; do diff $file test/$file; done
