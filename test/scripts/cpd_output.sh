#!/bin/bash

# currently existing amount of duplications in code
EXISTING_DUPLICATIONS=85

if [[ $(cat cpd.txt | wc -l) -gt 0 ]]; then
  # if the file contains any lines
  DUPLICATES=$(cat cpd.txt | grep "duplication in the following files" | wc -l)
  echo "----------------------------------------------------"
  echo "cpd found the following amount of duplications: " $DUPLICATES
  if [[ $DUPLICATES -gt $EXISTING_DUPLICATIONS ]]; then
    # fail if additional duplications were introduced
    echo "new duplications introduced, please fix"
    exit 1
  elif [[ $DUPLICATES -lt $EXISTING_DUPLICATIONS ]]; then
    # also fail if less duplications exist now, as a reminder to update the reference number
    echo "number of duplicates decreased, please update test/scripts/cpd_output.sh"
    exit 1
  else
    exit 0
  fi;
else
  exit 0
fi;
