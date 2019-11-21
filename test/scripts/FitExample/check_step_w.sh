 #!/bin/bash
diff -w -I "HistoAddress" -I "Opened input file" -I "Pruning" -I "libSM.so" -I "libASImage" -I "png file FitExample/Pruning.png" LOG_w test/logs/FitExample/LOG_w && diff FitExample/PruningText.txt test/reference/FitExample/PruningText.txt
