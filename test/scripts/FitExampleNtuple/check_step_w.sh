 #!/bin/bash
diff -w -I "HistoAddress" -I "Opened input file" -I "Pruning" -I "libSM.so" -I "libASImage" -I "png file FitExampleNtuple/Pruning.png" LOG_NTUPLE_w test/logs/LOG_NTUPLE_w && diff FitExampleNtuple/PruningText.txt test/FitExampleNtuple/PruningText.txt
