 #!/bin/bash
diff -w -I "HistoAddress" -I "Opened input file" -I "Pruning" -I "libSM.so" -I "libASImage" -I "png file FitExample/Pruning.png" LOG_w test/logs/LOG_w && diff FitExample/PruningText.txt test/FitExample/PruningText.txt
