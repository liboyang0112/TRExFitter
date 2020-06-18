#!/bin/bash
diff -w FitExampleUnfolding/Fits/FitExampleUnfolding.txt test/reference/FitExampleUnfolding/Fits/FitExampleUnfolding.txt && diff -w FitExampleUnfolding/Fits/UnfoldedResults.txt test/reference/FitExampleUnfolding/Fits/UnfoldedResults.txt
