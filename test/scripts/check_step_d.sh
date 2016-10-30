#!/bin/bash
diff LOG_d test/logs/LOG_d && for file in `ls FitExample/Tables/*.txt`; do diff $file test/$file; done
