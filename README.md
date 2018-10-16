# TRExFitter   [![build status](https://gitlab.cern.ch/TRExStats/TRExFitter/badges/master/build.svg "build status")](https://gitlab.cern.ch/TRExStats/TRExFitter/commits/master)

This package provides a framework to perform profile likelihood fits. In addition to that, many convenient features are available. TRExFitter was previously also known as TtHFitter. Here are a few important references to make use of:

* [TRExFitter twiki page](https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TtHFitter) for additional documentation and many references to further details
* [TRExFitter JIRA](https://its.cern.ch/jira/projects/TTHFITTER/summary>) (sign up to the mailing list in case you cannot access the JIRA)
* TRExFitter mailing list: [atlas-phys-stat-tthfitter](https://e-groups.cern.ch/e-groups/EgroupsSubscription.do?egroupName=atlas-phys-stat-tthfitter)


## Table of Contents
1.  [Getting the code](#getting-the-code)
2.  [Setup](#setup)
3.  [How to](#how-to)
4.  [Config File](#config-file)
    * [`Job` block options](#job-block-options)
    * [`Fit` block options](#fit-block-options)
    * [`Limit` block options](#limit-block-options)
    * [`Significance` block options](#significance-block-options)
    * [`Options` block options](#options-block-options)
    * [`Region` block options](#region-block-options)
    * [`Sample` block options](#sample-block-options)
    * [`ShapeFactor` block options](#shapefactor-block-options)
    * [`Systematic` block options](#systematic-block-options)
5.  [Command line options](#command-line-options)
6.  [Ranking Plot](#ranking-plot)
7.  [Grouped Impact](#grouped-impact)
8.  [Multi-Fit](#multi-fit)
    * [Multi-Fit `Job` block options](#multi-fit-job-block-options)
    * [Multi-Fit `Fit` block options](#multi-fit-fit-block-options)
9.  [Input File Merging with hupdate](#input-file-merging-with-hupdate)
10. [Output Directories Structure](#output-directories-structure)
11. [ShapeFactor example](#shapefactor-example)
12. [Replacement file](#replacement-file)
13. [TRExFitter package authors](#trexfitter-package-authors)



## Getting the code
To get the code, use the following command:
```
git clone ssh://git@gitlab.cern.ch:7999/TRExStats/TRExFitter.git
```
To get a specific tag, do the following:
```
cd TRExFitter && git checkout <tag number> && cd -
```



## Setup
To setup just use the script (from any location):
```
  source setup.sh
```
(should work on any machine with access to cvmfs - provided that nothing else is set-up previously)

To compile:

1. Create a new build directory with `mkdir build && cd build`
2. Run cmake with `cmake ../`
3. Compile the code with `cmake --build ./`
4. The binary file will appear in `bin/` directory

(this will take as main code the file `util/trex-fitter.C`)

The setup script also adds a path to the binary into your PATH and you can execute the code with `trex-fitter`

**IMPORTANT!** For the first time use you need to type `git submodule init` followed by `git submodule update`.
Every time the submodules change, you need to run `git submodule update`.

**Tip:** To recompile the code, directly from the main directory, a simple command is:
```
    cd build/ && cmake --build ./ && cd ../
```
or simply using the alias defined in the `setup.sh` script:
```
    trex-make
```


## How to
To run the code, after compiling (see [Setup](#setup)), use the command:
```
trex-fitter  <action(s)>  <config file>  [<options>]
```
The configuration file (`<config file>`) is a text file containing all the information on the definition of samples and fit regions, including all the fit and draw options.
By default, the file  `config/myFit.config`  is loaded.
See the section [Config File](#config-file) for more details.
Take a look at `config/myFit.config` or `config/ttH2015.config` to see some example config files.
Most of the time, the only file the user has to modify to obtain their desired results is the configuration file.

The only mandatory argument, `<action(s)>`, tells TRExFitter which operation(s) to perform.
The possible operations are defined in the main file (e.g. `util/trex-fitter.C`).
For instance, if you use the default file `util/trex-fitter.C`, the available actions are:

| **Option** | **Action** |
| ---------- | ---------- |
| `h` | read input histograms (valid only if the proper option is specified in the config file) |
| `n` | read input ntuples (valid only if the proper option is specified in the config file) |
| `w` | create the RooStats xmls and workspace |
| `f` | fit the workspace |
| `l` | calculate exclusion limit |
| `s` | calculate significance |
| `d` | draw pre-fit plots |
| `p` | draw post-fit plots |
| `a` | draw separation plots |
| `r` | draw ranking plot (see [Ranking Plot](#ranking-plot)) |
| `b` | re-run smoothing (in the future also rebinning) |
| `m` | multi-fit (see [Multi-Fit](#multi-fit)) |
| `i` | grouped impact evaluation (see [Grouped Impact](#grouped-impact)) |

New optional argument: `<options>`.
It is a string (so make sure to use " or ' to enclose the string if you use more than one option) defining a list of options, in the form:
```
"<option1>=<value1>,<value2>,...:<option2>=..."
```
See the section [Command line options](#command-line-options) below.



## Config File

The structure of the file should be the following:
```
<ObjectType>: <ObjectName>
  <ObjectProperty>: <Value>
  <ObjectProperty>: <Value>
  ...

<ObjectType>: <ObjectName>
  <ObjectProperty>: <Value>
  <ObjectProperty>: <Value>
  ...

...
```
NB: note the **blank** line between the objects!

The file should contain:
  * exactly one object of type `Job`
  * exactly one object of type `Fit`
  * exactly one object of type `Limit`
  * at least one object of type `Sample`
  * at least one object of type `Region`
  * any number of objects of type `Systematic` (even 0 is ok)
  * any number of objects of type `NormFactor` (even 0 is ok)

Note that each object should have unique `<ObjectName>`.

At the beginning of TRExFitter execution, the config file used will be checked against a reference file. The reference files for single and multi-fits are `jobSchema.config` and `multiFitSchema.config`, respectively. These files specify which options are allowed per block, and how the arguments should look like.

For each object type (also called "block"), here is the list of available properties:

### `Job` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| Label                        | the label which will be shown on plots |
| POI                          | the name of the parameter of interest; this should correspond to a NormFactor defined below |
| ReadFrom                     | can be HIST or NTUP; default is HIST |
| HistoPath(s)                 | valid only for option HIST above is selected; it's the path(s) where the input root files containing the histograms are stored |
| HistoFile(s)                 | valid only for option HIST; it's the file name(s) where the input root files containing the histograms are stored |
| HistoName(s)                 | valid only for option HIST; it's the histogram name(s) to read from the file(s) |
| NtuplePath(s)                | valid only for option NTUP; it's the path(s) where the input root files containing the ntuples are stored |
| NtupleFile(s)                | valid only for option NTUP; it's the file names(s) where the input root files containing the ntuples are stored |
| NtupleName(s)                | valid only for option HIST; it's the tree name(s) to read from the file(s) |
| MCweight                     | only for option NTUP; string defining the weight (for MC samples only) |
| Selection                    | only for option NTUP; string defining the selection |
| Lumi                         | value to scale all the "NormalizedByTheory" samples |
| LumiScale                    | additional value to scale 'after' histogram creation (for fast scaling) IMPORTANT: use it only if you know what you are doing!! |
| SystPruningShape             | Lower threshold to remove a shape systematic from the fit/limit (suppression is done per sample and per region) (e.g.: 0.02 for 2%) |
| SystPruningNorm              | Lower threshold to remove a normalisation systematic from the fit/limit (suppression is done per sample and per region) (e.g.: 0.02 for 2%) |
| SystLarge                    | all systematics above this threshold will be flagged in the pruning plot) (e.g. 0.4 will flag systematics that are larger than 40%) |
| IntCodeOverall               | interpolation code used for the normalization component of systematics (should match the one used in RooStats) |
| IntCodeShape                 | interpolation code used for the shape component of systematics (should match the one used in RooStats) |
| MCstatThreshold              | by default, the MC stat uncertainty is included in the fit (and to the plots); a NP will be added for each bin with an MC stat uncertainty > this threshold (relative) if the option is set to a float (default: no threshold); can also set to NONE in order to disable MC stat uncertainty completely |
| MCstatConstraint             | constraint used for MC stat uncertainties, can be set to 'GAUSSIAN' or 'POISSON' (default) |
| DebugLevel                   | 0 = prints only Warning and Errors, 1 = additionally prints Info messages, 2 = additionally prints Debug messages, >2 additionally prints Verbose messages. For option <2 RooFit/Roostats messages will be heavily suppressed |
| Logo                         | is set to TRUE will print the TRExFitter logo |
| PlotOptions                  | a set of options for plotting:<br>&nbsp; &nbsp; **YIELDS**: if set, the legend will be one-column and will include the yields; otherwise two-columns and no yields<br>&nbsp; &nbsp; **NORMSIG**: add normlised signal to plots<br>&nbsp; &nbsp; **NOSIG**: don't show signal in stack<br>&nbsp; &nbsp; **OVERSIG**: overlay signal (not normalised)<br>&nbsp; &nbsp; **CHI2**: the chi2/ndf and chi2 prob will be printed on each plot, provided that the option GetChi2 is set<br>&nbsp; &nbsp; **PREFITONPOSTFIT**: draw a dashed line on the postfit plot that indicates the sum of prefit background<br>&nbsp; &nbsp; **NOXERR**: removes the horizontal error bars on the data and the ratio plots |
| PlotOptionsSummary           | the same as PlotOptions but for the summary plot (if nothing is specified, PlotOptions is used) |
| TableOptions                 | a set of options for tables:<br>&nbsp; &nbsp; **STANDALONE**: default! If not set, no "\begin{document}"<br>&nbsp; &nbsp; **FOOTNOTESIZE**: -> \footnotesize <br>&nbsp; &nbsp; **LANDSCAPE**: -> \begin{landscape} |
| SystControlPlots             | if set to TRUE, plots showing the shape effect of a given systematic (before and after smoothing/symmetrisation) will be produced |
| SystDataPlots                | if set to TRUE, plots showing the shape effect of a given systematic (before and after smoothing/symmetrisation) on top of the nominal sum of samples will be produced. Data are then plotted in the ratio. If the option is set to "fillUpFrame", data will also be plotted in the upper frame. |
| CorrelationThreshold         | Threshold used to draw the correlation matrix (only systematics with at least one correlation larger than than draw) (0.05:5%) |
| SignalRegionsPlot            | list of regions to put in SignalRegionsPlot and PieChartPlots; use "EMPTY" to put an empty entry, "ENDL" to specify end of line. This specifies the order of regions plotted in signal region S/B plots and pie chart plots, as well as number of regions per row. |
| HistoChecks                  | NOCRASH: means that if an error is found in the input histograms, the code continues (with only warnings) -- default leads to a crash in case of problem |
| LumiLabel                    | label for luminosity to be put on plots |
| CmeLabel                     | label for center-of-mass energy to be put on plots |
| SplitHistoFiles              | set this to TRUE to have histogram files split by region (useful with many regions and/or run in parallel) |
| BlindingThreshold            | bins with S/B > this number will be blinded |
| KeepPrefitBlindedBins        | if set to TRUE, and if pre-fit an post-fit plots are produced together ("dp" option) pre-fit blinding is kept in post-fit plots |
| RankingMaxNP                 | max number of NP to show in ranking plot |
| RankingPlot                  | NP categories in gammas or systs, if set to Systs(Gammas) then plot only systs(Gammas) in ranking, default produce plot for systs+gammas, can also set to all to have the 3 plots. |
| ImageFormat                  | png, pdf or eps |
| StatOnly                     | the code ignores systematics and MC stat uncertainties from all computations (limits, significances, fit, ...); need to re-create ws in case of limit and significance |
| SystErrorBars                | TRUE by default to add stat error bars to syst variations in syst plots, set to FALSE to disable |
| SummaryPlotRegions           | list of regions to be shown in summary plot (useful to specify a custom order) |
| FixNPforStatOnly             | if set to TRUE, when running stat-only (with either of the two options) also the norm factors other than the POI are kept fixed |
| InputFolder                  | specify it to read fit input histograms from a different directory than `<jobName>/Histograms/` |
| InputName                    | specify it to read fit input histograms from files with different name than `<jobName>_blabla.root` |
| OutputDir                    | specify it to write everything in a different directory than `<jobName>` |
| WorkspaceFileName            | if specified, an external ws can be used as input for fitting (not 100% supported) |
| KeepPruning                  | if set to TRUE, the first time the ws is created (option w) a Pruning.root file is created under `<jobName>/` and used for future operations to skip pruned systematics (makes operations much faster in case many syst are pruned) |
| AtlasLabel                   | to specify Internal, Preliminary, etc... |
| CleanTables                  | if set to TRUE, a cleaned version of the tex tables is created (basically removing the "#") - to be expanded |
| SystCategoryTables           | if set to TRUE, additional syst tables with systematics grouped by category are created |
| SummaryPlotYmax              | if set, it will force the summary plot to use this value as max y-maxis value |
| SummaryPlotYmin              | if set, it will force the summary plot to use this value as min y-maxis value |
| RatioYmax                    | if set, it will specify the max of the range of the ratio plots |
| RatioYmin                    | if set, it will specify the min of the range of the ratio plots |
| RatioYmaxPostFit             | if set, it will specify the max of the range of the ratio plots, for post-fit only |
| RatioYminPostFit             | if set, it will specify the min of the range of the ratio plots, for post-fit only |
| CustomAsimov                 | if set, the workspace will be created with an AsimovData built according to Sample->`AsimovReplacementFor` option (see below) instead of data |
| RandomPOISeed                | if set to a >= 0 number, the signal sample(s) to which the POI is assigned get scaled by a random number generated starting from this seed, just before the ws creation; if the same seed is used in the cofig, post-fit plots will show consistent results (i.e. before post-fit drawing the POI is scaled by the same number) |
| GetChi2                      | if set to TRUE (or STAT+SYST), for pre- and post-fit plots the extended chi2 test is done, and results are printed on the screen for each plot when running d and/or p; can be set to STAT (or STAT-ONLY) for stat-only chi2 |
| TtresSmoothing               | if set to TRUE, the systematic uncertainty smoothing will use the ttbar resonances convention for the smoothing. The Smoothing parameter in the Systematics area can be set to 40 to treat the systematic uncertainty as correlated with the nominal or 400 to treat it as uncorrelated with the nominal. |
| SmoothingOption              | Choose which smoothing option to use, allowed parameters are: MAXVARIATION (default), TTBARRESONANCE, COMMONTOOLSMOOTHMONOTONIC, COMMONTOOLSMOOTHPARABOLIC, KERNELRATIOUNIFORM, KERNELDELTAGAUSS or KERNELRATIOGAUSS. |
| UseGammaPulls                | if set to TRUE, the fit results in terms of gamma parameter pulls, constraints and correlations are propagated to the post-fit plots, when possible (i.e. not for validation plots of course) |
| GuessMCStatEmptyBins         | For bins with negative yields, the yield is corrected to 1e-06. If a stat. uncertainty on that bin was defined, it is kept. If the stat. uncertainty was however 0, then this option kicks in. If it is set to TRUE (default), the smallest stat. uncertainty from any other bin in the distribution will be used. If set to FALSE, or if all other bins have no stat. uncertainty defined, a 100% uncertainty (of 1e-06) is applied instead. |
| CorrectNormForNegativeIntegral | By default, if there are samples with negative yields in some bins, the total integral over that sample will be re-scaled after the yield in the negative bins was corrected to 1e-6, such that the total yield across the sample in this region is consistent with the yield before fixing the negative bins. If this option is set to TRUE (by default, it is FALSE), then this re-scaling will also happen if the total yield initially was negative. This can lead to unexpected behavior (expert option, proceed with caution). |
| MergeUnderOverFlow           | if set to TRUE, the underflow content of each histogram is added to the first bin and the overflow to the last one (default is FALSE for HIST inputs and TRUE for NTUP inputs) |
| DoSummaryPlot                | if set to FALSE, no summary plot is created |
| DoMergedPlot                 | if set to TRUE, a merged plot of all the region groups specified with the RegionGroups option is created |
| DoTables                     | if set to FALSE, no tables are created |
| DoSignalRegionsPlot          | if set to FALSE, no signal regions plot is created |
| DoPieChartPlot               | if set to FALSE, no background composition pie-chart plot is created |
| CustomFunctions              | list of .C files with definition and implementation of functions to be used in strings defining selections or weights (see this link: https://wiki.physik.uzh.ch/lhcb/root:ttreedraw, notice that the file and function names should match and that all the arguments of the function should have default values) |
| SuppressNegativeBinWarnings  | If set to TRUE, will suppress warning messages about negative or 0 content in bins |
| Bootstrap                    | (only works with NTUP inputs) if set, the bootstrap method wil be used; the argument should be a string like `bsWeight(x,eventNumber,mcChannelNumber)`, where `bsWeight` should be loaded with `CustomFunctions: "bsWeight.C"` and eventNumber and mcChannelNumber should be existing branches for all the MC ntuples; then, to produce the i-th bootstrap pseudo-experiment, or to run on it (e.g. to perform a fit) the command-line option `BootstrapIdx=<i>` should be given, with `<i>=0,1,2,3...` |
| RunROOTMacros                | If set to True will run ROOT macros for limits and significa, otherwise (default) will run version which is compiled and has updated messaging. The functunality is the same. |
| DecorrSysts                  | comma-separated list of systematics which you want to decorrelate from another channel (this is don by automatically attaching a suffix to the NormFactor for each of them); can use wildcards |
| DecorrSuff                   | the suffix to attach when using DecorrSysts |
| RegionGroups                 | groups specified here will cause additional yield tables to be created per group, and also merged plots per group if DoMergedPlot is set to TRUE |
| ReplacementFile              | allows usage of placeholders in the config, which will be overwritten by values provided in an external file; see [Replacement file](#replacement-file) section |
| Suffix                       | added to file names of plots, workspace, fit results etc. (equivalent to command line option) |
| SaveSuffix                   | added to file name of histograms, for usage with hupdate (equivalent to command line option) |
| HideNP                       | comma-separated list of nuisance parameters to be excluded from pull plots and correlation matrix |
| SummaryPlotLabels            | DEPRECATED - labels to be used per region group in summary plot (only if FourTopStyle is set) |
| SummaryPlotValidationRegions | regions to be included in validation region summary plot (default: all) |
| SummaryPlotValidationLabels  | DEPRECATED - labels to be used per set of regions in validation-region summary plot (only if FourTopStyle is set) |
| SmoothMorphingTemplates      | if set to TRUE (default is FALSE), the templates used for morphig are forced to have linear dependence on the morphing parameter, bin-by-bin (plots are produced per bin, in the Morphing directory) |
| PropagateSystsForMorphing    | if set to TRUE (default is FALSE), the non-nominal templates inherit all the systematics from the nominal template; NB: the nominal template is determined by the nominal value of the norm factor in the config |
| SummaryPrefix                | adds a prefix to summary and merge plots |
| AllowWrongRegionSample       | Can be TRUE or FALSE (default). When set to TRUE code will print only warnings when chosen samples or regions for various options are not defined. When set to FALSE the code will print errors and stop when the samples/regions are not defined. |
| POIPrecision                 | Integer value N, N >=1 and N <=5. Will tell the code to use N decimal places for norm facotr mean value and uncertainty. Default is 2 |
| RankingPOIName               | Custom name for the POI for ranking plots. Default is `#mu` |
| UseGammasForCorr             | If set to `TRUE` will add gammas into correlation matrix plot. Default is `FALSE` |
| UseATLASRounding             | If set to `TRUE` will use PGD/ATLAS rounding to yield tables (both .txt and .tex) |
| UseATLASRoundingTxt          | If set to `TRUE` will use PGD/ATLAS rounding to yield tables (only .txt) |
| UseATLASRoundingTex          | If set to `TRUE` will use PGD/ATLAS rounding to yield tables (only .tex) |

### `Fit` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| FitType                      | can be SPLUSB (default) or BONLY to fit under the s+b or the b-only hypothesis |
| FitRegion                    | can be CRSR (default) or CRONLY to fit considering both signal and control regions in the fit, or only control regions. You can also specify a comma-separated list of regions to use in the fit |
| FitBlind                     | specify is real data or Asimov data should be used in the fit (TRUE or FALSE). By default, fit are NOT blind. |
| POIAsimov                    | value of the parameter of interest in the AsimovDataset used in the fit |
| NPValues                     | values of the nuisance parameters used to build the Asimov. Coma-separated list of NP:value (e.g. alpha_ttbarbb_XS:1,alpha_ttbarbcc_XS:1.5) |
| FixNPs                       | values of the nuisance parameters used to be fixed in the fit. Coma-separated list of NP:value (e.g. alpha_ttbarbb_XS:1,alpha_ttbarbcc_XS:1.5) |
| doLHscan                     | comma separated list of names of the POI or NP from which you want to produce the likelihood scan, if first element of the list is "all" then all systematics are profiled |
| LHscanMin                    | minimum value for the LH scan on x-axis (default it Norm min) |
| LHscanMax                    | maximum value for the LH scan on x-axis (default is Norm max) |
| LHscanSteps                  | number of steps on the LH scan (default is 30) |
| UseMinos                     | comma separated list of names of the POI and/or NP for which you want to calculate the MINOS errors, if first element of the list is "all" then the MINOS errors is calculated for all systematics and POIs |
| SetRandomInitialNPval        | useful to set this to >0 (e.g. 0.1) to help convergence of Asimov fits |
| SetRandomInitialNPvalSeed    | seed used to determine initial NP settings in minimization process if SetRandomInitialNPval option is enabled |
| NumCPU                       | specify the number of CPU to use for the minimization (default = 1) |
| StatOnlyFit                  | if specified, the fit will keep fixed all the NP to the latest fit result, and the fit results will be saved with the `_statOnly` suffix (also possible to use it from command line) |
| GetGoodnessOfFit             | set to TRUE to get it (based on chi2 probability from comparison of negative-log-likelihoods) |
| DoNonProfileFit              | [EXPERIMENTAL] if set to TRUE (default is FALSE), instead of the fit profilig the sysyetmatics, a set of stat-only fits will be performed, on an Asimov data-set created with one syst variation at a time |
| FitToys                      | [EXPERIMENTAL] if set to N > 0, N stat-ony toys are generated and fitted |
| ToysHistoMin                 | If FitToys is used, set minimum on the output toys histogram X axis |
| ToysHistoMax                 | If FitToys is used, set maximum on the output toys histogram X axis |
| ToysHistoNbins               | If FitToys is used, set number of bins for toys histogram output |
| TemplateInterpolationOption  | Option only for morping, tells the code which interpolation between the templates is used. Three possible options are available: LINEAR(default)/SMOOTHLINEAR/SQUAREROOT. All of these options basically use linear interpolation but SMOOTHLINEAR approximates it by integral of hyperbolic tangent and SQUAREROOT approximates it by $`\sqrt{x^2+\epsilon}`$ to achieve smooth transitions (first derivative) between the templates |

### `Limit` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| LimitType                    | can be ASYMPTOTIC or TOYS (the latter is not yet supported) |
| LimitBlind                   | can be TRUE or FALSE (TRUE means that ALL regions are blinded) |
| POIAsimov                    | value of the POI to inject in the Asimov dataset in LimitBlind is set to TRUE |
| SignalInjection              | if set to TRUE, expected signal with signal injection is evaluated |

### `Significance` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| SignificanceBlind            | can be TRUE or FALSE (TRUE means that ALL regions are blinded) |
| POIAsimov                    | value of the POI to inject in the Asimov dataset in SignificanceBlind is set to TRUE |

### `Options` block options:
  * additional options, accepting only float as arguments - useful for adding your functionalities & flags in a quick way, since they need minimal changes in the code) ...

### `Region` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| VariableTitle                | it's the label which will be displayed on the x-axis in the plots |
| Label                        | it's the label which will be showed on the plots and specifies which region is shown |
| TexLabel                     | label for tex files |
| ShortLabel                   | same as above, but a shorter version for plots with smaller available place |
| LumiLabel                    | label for luminosity to be put on plots |
| CmeLabel                     | label for center-of-mass energy to be put on plots |
| LogScale                     | set it to TRUE to have log-scale when plotting this region |
| HistoFile(s)                 | only for option HIST, the file name (or names, comma-separated) to be used |
| HistoName(s)                 | only for option HIST, the histogram name (or names, comma-separated) to be used |
| HistoPathSuff(s)             | only for option HIST, the path suffix (or suffixes, comma-separated) where to find the histogram files for this region |
| Variable                     | only for option NTUP, the variable (or expression) inside the ntuple to plot can define a variable as X|Y to do the correlation plot between X and Y |
| VariableForSample            | only for option NTUP, allows to set exceptions for Variable. This is a very useful feature when using TRF only in some samples. Comma-separated list of sample:variable (e.g. wjets:met_met/1e3,zjets:Mbbb/1e). |
| Selection                    | only for option NTUP, the selection done on the ntuple for this region |
| NtupleName(s)                | only for option NTUP, the name (or names, comma-separated) of the tree for this region |
| NtupleFile(s)                | only for option NTUP, the file (or files, comma-separated) of the tree for this region |
| NtuplePathSuff(s)            | only for option NTUP, the path sufix (or suffixes, comma-separated) where to find the ntuple files for this region |
| MCweight                     | only for option NTUP, the additional weight used in this region (for MC samples only) |
| Rebin                        | if specified, the histograms will be rebinned merging N bins together, where N is the argument (int) |
| Binning                      | if specified, the histograms will be binned according to the new binning specified, in the form like (0,10,20,50,100). If option AutoBin is set, use algorithms/functions or define the binning. Example - Binning: "AutoBin","TransfoD",5.,6. (TransfoF also available, 5. and 6. are parameters of the transformation). If used in background region and zSig!=0 (first parameter, =0 gives flat background) then need a comma separated list of backgrounds to use instead of signal to compute the binning. |
| Rebinning                    | if specified, the histograms will be rebinned according to the new binning specified, in the form like (0,10,20,50,100). Differently from the BInning option, this one performs the rebinning aftre the orginal histograms are created. This means that this option can changed (or removed) before running the b step. |
| BinWidth                     | if specified, two things are done: this number is used to decorate the y axis label and the bin content is scaled for bins with a bin width different from this number |
| BinLabels                    | if specified, bin labels are set according to provided comma separated list (list length must be equal to number of bins) |
| Type                         | can be SIGNAL, CONTROL or VALIDATION; used depending on Fit->FitType; if VALIDATION is set, the region is never fitted; default is SIGNAL |
| DataType                     | ASIMOV or DATA. Is Asimov is set, the limits and significances are computed without taking into account the data in these region, but a projection of the fit performed in the regions with DATA |
| Ymax                         | if set, it will force the plot to use this value as max y-maxis value |
| Ymin                         | if set, it will force the plot to use this value as min y-maxis value |
| RatioYmax                    | if set, it will specify the max of the range of the ratio plot for this region only |
| RatioYmin                    | if set, it will specify the min of the range of the ratio plot for this region only |
| RatioYmaxPostFit             | if set, it will specify the max of the range of the ratio plot for this region only, for post-fit only |
| RatioYminPostFit             | if set, it will specify the min of the range of the ratio plot for this region only, for post-fit only |
| DropBins                     | allows to specify a comma-separated list of bins to set to 0 (both for data and prediction), starting from 0 for the index |
| Group                        | if specified, regions of the same group appear together in several places, see RegionGroups option |
| YaxisTitle                   | title of y-axis used for plots of the region |
| YmaxScale                    | scales range of y-axis (default: 2.0, meaning the maximum axis value is twice the largest yield in any bin) |
| Ymax                         | maximum value on y-axis |
| SkipSmoothing                | if smoothing of nominal samples is used, this option can be used to disable smoothing per region (default: FALSE) |

### `Sample` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| Type                         | can be SIGNAL, BACKGROUND, DATA or GHOST; default is BACKGROUND; GHOST means: no syst, not drawn, not propagated to workspace |
| Title                        | title shown on the legends |
| TexTitle                     | title shown on tex tables |
| Group                        | if specified, sample will be grouped with other samples with same group and this label will be used in plots |
| HistoName(s)                 | valid only for option HIST; name(s) of histogram to read |
| HistoFile(s)                 | valid only for option HIST; which root file(s) to read (excluding the suffix ".root"); this will be combined with Fit->HistoPath to build the full path |
| HistoPath(s)                 | valid only for option HIST; it's the path(s) where the input root files containing the histograms are stored |
| HistoNameSuff(s)             | valid only for option HIST; suffix(es) for the name of histogram(s) to read |
| HistoFileSuff(s)             | valid only for option HIST; suffix(es) for the file name of histogram(s) to read |
| HistoPathSuff(s)             | valid only for option HIST; suffix(es) for the path of histogram(s) to read |
| NtupleName(s)                | valid only for option NTUP; name(s) of tree to read |
| NtupleFile(s)                | valid only for option NTUP; it's the file name(s) where the input ntuples are stored |
| NtuplePath(s)                | valid only for option NTUP; it's the path(s) where the input root files containing the ntuples are stored |
| NtupleNameSuff(s)            | valid only for option NTUP; suffix(es) for the name of tree to read |
| NtupleFileSuff(s)            | valid only for option NTUP; suffix(es) for the file name(s) of tree to read |
| NtuplePathSuff(s)            | valid only for option NTUP; suffix(es) for the path(s) of tree to read |
| FillColor                    | histogram fill color (not valid for data) |
| LineColor                    | histogram line color |
| NormFactor                   | NormalisationFactor (free parameter in the fit); in the format \<name\>,nominal,min,max |
| ShapeFactor                  | ShapeFactor added |
| NormalizedByTheory           | set it to FALSE for data-driven backgrounds (MCweight, Lumi and LumiScale from Job and Region will be ignored) |
| MCweight                     | only for option NTUP, the additional weight used in this sample (for all types of samples!! Not only MC) |
| Selection                    | valid only for option NTUP; additional selection for this region |
| Regions                      | set this to have the sample only in some regions |
| Exclude                      | set this to exclude the sample in some regions |
| LumiScale(s)                 | set this to scale the sample by a number; if more numbers are set, use a different one for each file / name / path... |
| IgnoreSelection              | if set, selection from Job and Region will be ignored |
| UseMCstat                    | if set to FALSE, makes the fitter ignore the stat uncertainty for this sample |
| UseSystematics               | has to be set to TRUE to allow systematics on the GHOST samples |
| MultiplyBy                   | if specified, each sample hist is multiplied bin-by-bin by another sample hist, in each of the regions |
| DivideBy                     | if specified, each sample hist is divided bin-by-bin by another sample hist, in each of the regions |
| AddSample(s)                 | if specified, each sample hist gets added bin-by-bin another sample hist, in each of the regions |
| SubtractSample(s)            | if specified, each sample hist gets subtracted bin-by-bin another sample hist, in each of the regions |
| Smooth                       | if set to TRUE, the nominal histograms are smoothed (based on TH1::Smooth but taking into account the original stat uncertainty) |
| AsimovReplacementFor         | only for GHOST samples; if set, the creation of custom Asimov data-set(s) is triggered; use as 'AsimovReplacementFor: "dataset","sample"', where "dataset" is the name of a custom Asimov dataset one wants to create (the same name will have to be set under Job->CustomAsimov in order to use it) and "sample" is the sample this GHOST sample will supersede |
| SeparateGammas               | if set to TRUE, the sample will not contribute to the overall gamma factors for MC stat, but a separate set of them will be added for this sample (through the SHAPE systematic technology); NB: you need to re-run at least the "b" step if you want to decorrelate the gammas on existing inputs (wf is not enough) |
| CorrelateGammasInRegions     | to be used only together with SeparateGammas; can be used to correlate MC stat across regions; example: "SR1:SR2,CR1:CR2:CR3" will use the same NP for the MC stat in each bin of SR1 and SR2 and in each bin of CR1, CR2 and CR3 |
| CorrelateGammasWithSample    | to be used only together with SeparateGammas; can be used to correlate MC stat of this sample with those of another sample (example usecase: when one sample is derived from another one through reweighting, and they are both used in the fit) |
| Morphing                     | add this to each template you have, to do a template fit / morphing; syntax is `<name-of-parameter>,<value-corresponding-to-this-template>`; the POI should be set to `<name-of-parameter>` |
| BuildPullTable               | if set to TRUE or NORM-ONLY, create tables showing the post-fit acceptance effect of nuisance parameter pulls for this sample, set to NORM+SHAPE to include the bin-by-bin effect |

### `NormFactor` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| Samples                      | comma-separated list of samples on which to apply the norm factor |
| Regions                      | comma-separated list of regions where to apply the norm factor |
| Exclude                      | comma-separated list of samples/regions to exclude |
| Title                        | title of the norm factor |
| Nominal                      | nominal value |
| Min                          | min value |
| Max                          | max value |
| Constant                     | set to TRUE to have a fixed norm factor |
| Category                     | major category to which the NormFactor belongs (instrumental, theory, ttbar, ...) |
| SubCategory                  | minor category for the NormFactor, used to evaluate impact on POI per SubCategory in "i" step, defaults to "NormFactors", do not use "Gammas", "FullSyst", or "combine" as SubCategory names (reserved for special functionality) |
| Expression                   | a way to correlate this norm factor with other norm factors (using AddPreprocessFunction); two argments, in the form `<expression>,<dependency>`, where `<dependency>` should contain the name(s) of the norm factor the expression depends on [example: `"1.-SigXsecOverSM","SigXsecOverSM"`] |

### `ShapeFactor` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| Samples                      | comma-separated list of samples on which to apply the shape factor |
| Regions                      | comma-separated list of regions where to apply the shape factor |
| Title                        | title of the shape factor |

### `Systematic` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| Samples                      | comma-separated list of samples on which to apply the systematic |
| Regions                      | comma-separated list of regions where to apply the systematic |
| Exclude                      | comma-separated list of samples/regions to exclude |
| ExcludeRegionSample          | comma-separated list of region:sample to exclude |
| Type                         | can be HISTO, OVERALL, SHAPE (this refers to the HistFactory Shape Systematic, i.e. uncorrelated bin-by-bin) or STAT (this refers to auto-creation of one systematic from stat uncertainty for each bin of corresponding region - DEPRECATED) |
| Title                        | title of the systematic (will be shown in plots) |
| StoredName                   | if specified, will be used to read and write histograms in the root files under Histograms/ intead of the syst name; useful to decorrelate without re-creating histograms |
| NuisancaParameter            | if specified, this will be given to RooStats instead of the syst name; useful (and recommended) way to correlate systematics |
| IsFreeParameter              | if set to TRUE, the constraint will be a flat one instead of Gaussian (use with caution) |
| Category                     | major category to which the systematic belongs (instrumental, theory, ttbar, ...): used to split pulls plot for same category |
| SubCategory                  | minor category for the systematic, used to evaluate impact on POI per SubCategory in "i" step, defaults to Category setting if it is used, otherwise defaults to "Uncategorised", do not use "Gammas", "FullSyst", or "combine" as SubCategory names (reserved for special functionality) |
| HistoPathUp                  | only for option HIST, for HISTO or SHAPE systematic: histogram file path for systematic up variation |
| HistoPathDown                | only for option HIST, for HISTO or SHAPE systematic: histogram file path for systematic down variation |
| HistoPathSufUp               | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic up variation |
| HistoPathSufDown             | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic down variation |
| HistoFileUp                  | only for option HIST, for HISTO or SHAPE systematic: histogram file name for systematic up variation |
| HistoFileDown                | only for option HIST, for HISTO or SHAPE systematic: histogram file name for systematic down variation |
| HistoFileSufUp               | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic up variation |
| HistoFileSufDown             | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic down variation |
| HistoNameUp                  | only for option HIST, for HISTO or SHAPE systematic: histogram name for systematic up variation |
| HistoNameDown                | only for option HIST, for HISTO or SHAPE systematic: histogram name for systematic down variation |
| HistoNameSufUp               | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram names for systematic up variation |
| HistoNameSufDown             | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram names for systematic down variation |
| NtuplePath(s)Up              | only for option NTUP, for HISTO or SHAPE systematic: ntuple file path(s) for systematic up variation |
| NtuplePath(s)Down            | only for option NTUP, for HISTO or SHAPE systematic: ntuple file path(s) for systematic down variation |
| NtuplePathSufUp              | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file paths for systematic up variation |
| NtuplePathSufDown            | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file paths for systematic down variation |
| NtupleFile(s)Up              | only for option NTUP, for HISTO or SHAPE systematic: ntuple file name(s) for systematic up variation |
| NtupleFile(s)Down            | only for option NTUP, for HISTO or SHAPE systematic: ntuple file name(s) for systematic down variation |
| NtupleFileSufUp              | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file names for systematic up variation |
| NtupleFileSufDown            | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file names for systematic down variation |
| NtupleName(s)Up              | only for option NTUP, for HISTO or SHAPE systematic: ntuple name(s) for systematic up variation |
| NtupleName(s)Down            | only for option NTUP, for HISTO or SHAPE systematic: ntuple name(s) for systematic down variation |
| NtupleNameSufUp              | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple names for systematic up variation |
| NtupleNameSufDown            | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple names for systematic down variation |
| SampleUp                     | if set, the syst variation will be built comparing the sample with another sample after all corrections are done; NB: can be used only if the syst affects one sample only |
| SampleDown                   | if set, the syst variation will be built comparing the sample with another sample after all corrections are done; NB: can be used only if the syst affects one sample only |
| WeightUp                     | only for option NTUP, for HISTO or SHAPE systematic: weight for systematic up variation; it will be multiplied with the `MCweight` defined in the Job and Region blocks, however it will not be multiplied with the `MCweight` defined in the sample block of the reference sample (for this, use `WeightSufUp` instead) |
| WeightDown                   | only for option NTUP, for HISTO or SHAPE systematic: weight for systematic down variation; it will be multiplied with the `MCweight` defined in the Job and Region blocks, however it will not be multiplied with the `MCweight` defined in the sample block of the reference sample (for this, use `WeightSufDown` instead) |
| WeightSufUp                  | only for option NTUP, for HISTO or SHAPE systematic: additional weight for systematic up variation (also multiplied with the `MCweight` acting on the nominal sample) |
| WeightSufDown                | only for option NTUP, for HISTO or SHAPE systematic: additional weight for systematic down variation (also multiplied with the `MCweight` acting on the nominal sample) |
| IgnoreWeight                 | only for option NTUP: if set, the corresponding weight (present in Job, Sample or Region) will be ignored for this systematic |
| Symmetrisation               | can be ONESIDED or TWOSIDED (...); for no symmetrisation, skip the line |
| Smoothing                    | smoothing code to apply; use 40 for default smoothing; for no smoothing, skip the line |
| OverallUp                    | for OVERALL systematic: the relative "up" shift (0.1 means +10%) |
| OverallDown                  | for OVERALL systematic: the relative "down" shift (-0.1 means -10%) |
| ScaleUp                      | for OVERALL, HISTO or SHAPE systematic: scale difference between "up" and nominal by a factor, or different factors for different regions (with the syntax `region1:1.2,region2:0.9`) |
| ScaleDown                    | for OVERALL, HISTO or SHAPE systematic: scale difference between "down" and nominal by a factor, or different factors for different regions (with the syntax `region1:1.2,region2:0.9`) |
| ReferenceSample              | if this is specified, the syst variation is evaluated w.r.t. this reference sample (often a GHOST sample) instead of the nominal, and then the relative difference is propagated to nominal; NOTE: also the overall relative difference is propagated |
| DropShapeIn                  | specify regions where you want the smoothing / pruning to be forced to drop the shape and keep only norm |
| DropNorm                     | the same as the previous one, but to drop the norm and keep only the shape |
| KeepNormForSamples           | list of samples (or sum of samples, in the form smp1+smp2), comma separated, for which the systematic gets shape only in each region |
| PreSmoothing                 | if set to TRUE, a TH1::Smooth-based smoothing is applied, prior to the usual smoothing (if set) |
| SubtractRefSampleVar         | if set to TRUE, the relative variation of the ReferenceSample will be linearly subtracted from the relative variation of each affected sample, for the same systematic - this is relevant e.g. for Full JER SmearingModel, where data would be the reference sample |
| HistoPathUpRefSample         | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file path for systematic up variation |
| HistoPathDownRefSample       | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file path for systematic down variation |
| HistoPathSufUpRefSample      | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram file names for systematic up variation |
| HistoPathSufDownRefSample    | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram file names for systematic down variation |
| HistoFileUpRefSample         | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file name for systematic up variation |
| HistoFileDownRefSample       | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file name for systematic down variation |
| HistoFileSufUpRefSample      | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram file names for systematic up variation |
| HistoFileSufDownRefSample    | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram file names for systematic down variation |
| HistoNameUpRefSample         | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram name for systematic up variation |
| HistoNameDownRefSample       | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram name for systematic down variation |
| HistoNameSufUpRefSample      | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram names for systematic up variation |
| HistoNameSufDownRefSample    | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram names for systematic down variation |
| NtuplePath(s)UpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple file path(s) for systematic up variation |
| NtuplePath(s)DownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple file path(s) for systematic down variation |
| NtuplePathSufUpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple file paths for systematic up variation |
| NtuplePathSufDownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple file paths for systematic down variation |
| NtupleFile(s)UpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple file name(s) for systematic up variation |
| NtupleFile(s)DownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple file name(s) for systematic down variation |
| NtupleFileSufUpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple file names for systematic up variation |
| NtupleFileSufDownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple file names for systematic down variation |
| NtupleName(s)UpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple name(s) for systematic up variation |
| NtupleName(s)DownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple name(s) for systematic down variation |
| NtupleNameSufUpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple names for systematic up variation |
| NtupleNameSufDownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple names for systematic down variation |
| Decorrelate                  | decorrelate systematic, can take values REGION (decorrelate across regions), SAMPLE (decorrelate across samples), SHAPEACC (decorrelate shape and acceptance effects) |



## Command line options
Currently the supported options are:

| **Option** | **Effect** |
| ---------- | ---------- |
| **Regions**       | to limit the regions to use to the list specified |
| **Samples**       | to limit the samples to use to the list specified |
| **Systematics**   | to limit the systematics to use to the list specified |
| **Signal**        | in case more than one SIGNAL sample is specified in your config file, you can specify which one you want to run on (for plots, workspace creation and fits/limits/significance) |
| **Exclude**       | to exclude certain Regions / Samples / Systematics |
| **Suffix**        | used for: plots, workspace, fit results, etc |
| **SaveSuffix**    | used for: saving histograms with a suffix (to be merged / renamed later, see [Input File Merging with hupdate](#input-file-merging-with-hupdate) section |
| **Update**        | if TRUE, the output .root file is updated, otherwise is overwrote |
| **StatOnlyFit**   | if TRUE, the same as Fit->StatOnlyFit |
| **StatOnly**      | if TRUE, no systematics nor MC stat uncertainties will be considered (equivalent to set StatOnly: TRUE in the Job block of the config), use `Systematics=NONE` instead to keep MC stat uncertainties |
| **Ranking**       | see [Ranking Plot](#ranking-plot) section |
| **FitResults**    | the specified fit results file will be used, for instance for post-fit plots (instead of the file `jobName/Fits/jobName.txt`) |
| **FitType**       | can be set to SPLUSB or BONLY to replace the option in the config file |
| **LumiScale**     | as the options in config file |
| **BootstrapIdx**  | see description of Bootstrap option in config (under Job) |
| **GroupedImpact** | see [Grouped Impact](#grouped-impact) section |

Note: the wild-card `*` is supported, but only as last character.
Example:
```
trex-fitter  n  config/ttH2015.config 'Regions=HThad_ge6jge4b:Exclude=BTag_*'
```



## Ranking Plot
* The ranking plot can be created in one go, with just the command line argument `r` (after having run the nominal fit `f`).
* Since this can take too much time (and memory), for complicated fits it's better to run it in several steps:
   by specifying the command-line option `Ranking=<name/index>`, one can produce the txt input for the ranking only for a specific line of the ranking, i.e. for a single NP (specified either through its name or index). Once all the needed txt files are created (e.g. in parallel through batch jobs) with the option `Ranking=plot` they are merged to create the final plot.

* Examples:
```
# this runs the ranking in one go
trex-fitter  r  <config>
#these commands will first create the inputs for the ranking one by one and then merge them in the plot
trex-fitter  r  <config> Ranking=Lumi
trex-fitter  r  <config> Ranking=JES1
trex-fitter  r  <config> Ranking=ttXsec
trex-fitter  r  <config> Ranking=plot
```



## Grouped Impact
* The command line argument `i` is used to evaluate the combined impact of groups of nuisance parameters on the POI.
* Specify groups using the `SubCategory` option (for Systematics and NormFactors).
* Two groups are defined by default: "Gammas" (MC stat. impact) and "FullSyst" (full systematics impact with statistical component subtracted).
* The impact is calculated by performing a fit where the nuisance parameters in the group are fixed to their best-fit values, and then the subtracting the resulting uncertainty on the POI in quadrature from the uncertainty from the nominal fit.
* The command line parameter `GroupedImpact` can be used to parallelize the impact calculations. If it is not specified, all existing groups are evaluated sequentially.
* The results are saved in `Fits/GroupedImpact*`.

Examples:

```
# evaluate impact of all groups sequentially
trex-fitter i <config>

# evaluate only the impact of Gammas
trex-fitter i <config> GroupedImpact="Gammas"
```

If the calculations are parallelized, combine the results by running the following at the end:

```
trex-fitter i <config> GroupedImpact="combine"
```



## Multi-Fit
The Multi-Fit functionality can be used to compare fit results or even to combine fit inputs from different configuration files / Jobs.

To use it you need a dedicated config file, with a structure similar to the usual ones. Example:
```
MultiFit: "myTopWS_multifit"
  Label: "My Label"
  Combine: FALSE
  Compare: TRUE
  CmeLabel: "13 TeV"
  LumiLabel: "85 pb^{-1}"
  ComparePOI: TRUE
  ComparePulls: TRUE
  CompareLimits: TRUE
  POIName: "SigXsecOverSM"
  POIRange: -10,30
  DataName: "obsData"
  CombineChByCh: TRUE

Fit: "CR"
  ConfigFile: config/myTopWS_CR.config
  Label: "CR-only"

Fit: "SR"
  ConfigFile: config/myTopWS_SR.config
  Label: "SR"
```
This config file can be run with the command line:
```
trex-fitter  m  config/myTopWS_multifit.config
```
  This will compare the fit results in terms of fitted NP, fitted POI and limits from the two config files specified. Notice that the fit and limits results have to be already available (they are not produced on the fly when running his multi-fit option).

To make a real combination, one needs to use the usual command options `w`, `f` and `l` together with the flag "Combine: TRUE" in the config above. Example:
```
trex-fitter  mwf  config/myTopWS_multifit.config
```
This will create a combined ws starting from the individual ws for the different regions in the two config files, and fit it.


### Multi-Fit `Job` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| Label            | the label which will be shown on plots |
| OutputDir        | the name of the output directory |
| LumiLabel        | the luminosity label that will be shown on the pltos |
| CmeLabel         | the center of mass energy label that will be shown on the plots |
| SaveSuf          | added to file name of histograms, for usage with hupdate (equivalent to command line option) |
| ShowObserved     | can be TRUE or FALSE, flag to turn on/off the observed values on the plots |
| LimitTitle       | the title for limit that will be shwon on the plots |
| POITitle         | the title of the POI that will be shown on X axis  |
| CompareLimits    | can be TRUE or FALSE, flag to compare to Limit values |
| ComparePOI       | can be TRUE or FALSE, flag to compare to POI values |
| ComparePulls     | can be TRUE or FALSE, flag to compare to pulls values |
| PlotCombCorrMatrix | can be set to TRUE or FALSE, flag to build correlation matrix from the combined systematics |
| Combine          | can be TRUE or FALSE, set to TRUE if you want to perfom actual combination (followed by `mwf`) |
| Compare          | can be TRUE or FALSE, set to TRUE if you want to compare values |
| StatOnly         | can be TRUE or FALSE, set to TRUE if the fits are stat only fits |
| IncludeStatOnly  | can be TRUE or FALSE, set to TRUE if you want to include stat only fits |
| POIName          | the name of the POI in the configs |
| POIRange         | the range of the chosen POI |
| LimitMax         | set maximum value for the limit |
| POIVal           | the value of the POI (for ASIMOV) |
| POIPrecision     | string, set precision of the POI |
| DataName         | can be "obsData", "asimovData", or custom string, if nothing is specified the observed data will be used |
| FitType          | can be SPLUSB or BONLY |
| SignalInjection  | can be TRUE or FALSE |
| CombineChByCh    | can be TRUE (default) or FALSE, set to TRUE to combine channel by channel |
| NPCategories     | comma separated list of NP categories |
| SetRandomInitialNPval | provide a float  |
| SetRandomInitialNPvalSeed | provide an int |
| NumCPU           | a number of CPU cores used for the fit |
| FastFit          | can be TRUE or FALSE |
| FastFitForRanking| can be TRUE or FALSE |
| NuisParListFile  |
| PlotSoverB       | if set to TRUE will plot signal over background plots |
| SignalTitle      | a title of the signal for the plots |
| FitResultsFile   | a name of the file with fit results |
| LimitsFile       | a name of the file with limits results |
| BonlySuffix      | a suffix of the background only fits |
| ShowSystForPOI   | can be TRUE or FALSE, set to TRUE if you want to show systematics for POI |
| GetGoodnessOfFit | can be TRUE or FALSE, set to TRUE to get chi2/NDF for the fit |
| doLHscan         | comma separeted list of NP(or POIs) to run LH scan, if first parameter is "all" it will be run for all NP |
| LHscanMin        | minimum value for the LH scan on x-axis (default it Norm min) |
| LHscanMax        | maximum value for the LH scan on x-axis (default is Norm max) |
| LHscanSteps      | number of steps on the LH scan (default is 30) |
| PlotOptions      | same as for "standard" fits |
| Logo             | can be TRUE or FALSE, use TRUE to show TRExFitter logo |
| DebugLevel       | set level of debug output |
| RunROOTMacros    | can be TRUE or FALSE, set to TRUE to run the common scripts in root interpreter in stead of running the directly compiled version (FALSE, default) |
| POILabel         | name of the POI shwon on plots, default is `#\mu` |
| POINominal       | value of the nominal (SM) prediction for POI, default is `1` |
| ShowTotalOnly    | If set to TRUE will show only total uncertainty on the POI plots. Default is FALSE |

### Multi-Fit `Fit` block options:

| **Option** | **Function** |
| ---------- | ------------ |
| Options          | additional options, accepting only float as arguments - useful for adding your functionalities & flags in a quick way, since they need minimal changes in the code) ... |
| Label            | the label of the values from this config that will be shown on the plots |
| LoadSuf          |
| ConfigFile       | the path to the config file that you want to combine/compare |
| Workspace        | the path to the worskapce that you want to combine/compare |
| ShowObserved     | can be TRUE or FALSE, set to TRUE to show the observed values of POI |
| FitResultsFile   | the path to the file with fit results |
| LimitsFile       | the path to the file with limits results |
| POIName          | the name of the POI |
| Directory        | the path to the directory |
| InputName        | the name of the input |



## Input File Merging with hupdate
A macro `hupdate` is included, which mimics hadd functionality, but without adding histograms if they have the same name.
This is useful for running different systematics in different steps (like different batch jobs) and then merging results afterwards.
`hupdate` is compiled automatically when using cmake. To explicitly request compilation, execute the following in the build folder:
```
make hupdate.exe
```
Example usage, combined with the usage of SaveSuffix:
```
make hupdate.exe
trex-fitter n ../config/ttH2015.config Systematics=BTag_B_NP1:SaveSuffix=_BTag_B_NP1
./build/bin/myFit.exe n ../config/ttH2015.config Exclude=BTag_B_NP1:SaveSuffix=_rest
./build/bin/hupdate.exe ../ttH2015/Histograms/ttH2015_HThad_4j2b_histos.root ttH2015/Histograms/ttH2015_HThad_4j2b_histos_rest.root ttH2015/Histograms/ttH2015_HThad_4j2b_histos_BTag_B_NP1.root
./build/bin/hupdate.exe ../ttH2015/Histograms/ttH2015_HThad_5j3b_histos_NEW.root ttH2015/Histograms/ttH2015_HThad_5j3b_histos.root ttH2015/Histograms/ttH2015_HThad_5j3b_histos_BTag_B_NP1.root
./build/bin/hupdate.exe ../ttH2015/Histograms/ttH2015_HThad_ge6jge4b_histos_NEW.root ttH2015/Histograms/ttH2015_HThad_ge6jge4b_histos.root ttH2015/Histograms/ttH2015_HThad_ge6jge4b_histos_BTag_B_NP1.root
trex-fitter dwf ../config/ttH2015.config
```



## Output Directories Structure
For each TRExFit object, a directory is created, with the same name as the Fit Name.
Inside this directory, at every step, some outputs are created, following the structure described above:

| **Folder** | **Content** |
| ---------- | ----------- |
| `Plots/`              | data/MC plots, pre- and post-fit, for all the Signal, Control and Validation regions, including the summary plots |
| `Tables/`             | tables in txt and tex format |
| `RooStats/`           | workspace(s) and the xmls |
| `Fits/`               | output from fits |
| `Limits/`             | outputs from the limit-setting code |
| `Significance/`       | outputs from the significance code |
| `Systematics/`        | plots for the syst variations |
| `Toys/`               | plots and ROOT files with pseudoexperiments output |
| `Histograms/`         | root file(s) with all the inputs |
| `LHoodPlots/`         | likelihood scan with respect to the specified parameter |



## ShapeFactor example
* The following scripts create example histograms in `exampleDataDriven` directory and execute `trex-fitter` using `config/dataDriven.config`
* The example contains a control region and signal region with two bins. The shape of one of the background samples is estimated using the ShapeFactor:
```
python makeDataDriven.py
python runDataDrivenExample.py
```
The results are in `JobDataDriven`



## Replacement file
You can define placeholders in your config file, which are replaced with values specified in an external file, which is read at the beginning of TRExFitter execution. This requires adding an additional option into your config, as part of the Job block:
```
ReplacementFile: path/to/file.txt
```
The replacement file should have the following structure:
```
# comment
XXX_placeholder: 0.1
XXX_another_placeholder: 0.2
% also a comment
```
Note that all placeholders must start with ``XXX``. In your config file, you can then refer to the placeholders like this:
```
Sample: "ttbar"
  MCweight: XXX_placeholder
```
If you would like to ensure that the replacement works correctly, set your `DebugLevel` to a minimum value of 1 and check the output of the framework.



## TRExFitter package authors
Managers:

* Michele Pinamonti [michele.pinamonti@gmail.com](michele.pinamonti@gmail.com)
* Loic Valery [loic.valery@cern.ch](loic.valery@cern.ch)

Development and support team:

* Alexander Held [alexander.held@cern.ch](alexander.held@cern.ch)
* Tomas Dado [tomas.dado@cern.ch](tomas.dado@cern.ch)
