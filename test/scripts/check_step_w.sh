#!/bin/bash
diff -I "HistoAddress" -I "Opened input filo" -I "Pruning" LOG_w test/logs/LOG_w && diff FitExample/PruningText.txt test/FitExample/PruningText.txt
