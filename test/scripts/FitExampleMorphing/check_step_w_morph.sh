 #!/bin/bash
diff -w -I "HistoAddress" -I "Opened input file" -I "Pruning" -I "libSM.so" -I "libASImage" -I "png file FitExampleMorphing/Pruning.png" LOG_MORPH_w test/logs/FitExampleMorphing/LOG_MORPH_w && diff FitExampleMorphing/PruningText.txt test/reference/FitExampleMorphing/PruningText.txt
