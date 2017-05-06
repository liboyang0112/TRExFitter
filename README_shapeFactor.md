# README shapeFactor#

### Adding ShapeFactor to TRexFitter ###

See [JIRA ticket TTHFITTER-115](https://its.cern.ch/jira/browse/TTHFITTER-115) for motivation

### Set up ###

Checkout the development  branch

```
git clone ssh://git@gitlab.cern.ch:7999/TRExStats/TRExFitter.git
cd TRExFitter
git checkout shapeFactor
source setup.sh
time make
```

### How to run it ###
Creates example histograms in `exampleDataDriven` directory and executes `myFit.exe` using `config/dataDriven.config`

```
python makeDataDriven.py
python runDataDrivenExample.py
```
The results are in `JobDataDriven`

To browse the differences of the branches per file one can do e.g. the following:

```
git config --global diff.tool meld
git difftool master shapeFactor

```

### Todo ###

* Check post-fit plots and tables
* Create separate NP plot for ShapeFactor