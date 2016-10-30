#!/bin/bash
echo "=> Checking the w step by comparing the logfiles compared to expectations"
diff -I "HistoAddress" -I "Opened input filo" -I "Pruning" LOG_w test/logs/LOG_w && diff FitExample/PruningText.txt test/FitExample/PruningText.txt && echo "toto"
echo "=> Done checking the w step"
