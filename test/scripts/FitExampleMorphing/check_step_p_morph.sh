#!/bin/bash
diff -w -I "libSM.so" -I "libASImage" -I "png file FitExampleMorphing" LOG_MORPH_p test/logs/FitExampleMorphing/LOG_MORPH_p && for file in `ls FitExampleMorphing/Tables/*postFit.txt`; do diff -w $file test/reference/$file; done
