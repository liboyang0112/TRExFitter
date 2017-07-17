Build status
---------
.. image:: https://gitlab.cern.ch/TRExStats/TRExFitter/badges/master/build.svg 
   :target: https://gitlab.cern.ch/TRExStats/TRExFitter/commits/master


Getting the code
---------
To get the code, use the following command::

  git clone ssh://git@gitlab.cern.ch:7999/TRExStats/TRExFitter.git

To get a specific tag, do the following::

  cd TRExFitter && git checkout <tag number> && cd -

Setup
---------
To setup just use the script::

  source setup.sh

(should work on any machine with access to afs - provided that nothing else is set-up previously)

To compile::

  make

(this will take as main code the file util/myFit.C)

or specify a custom main file as argument::

  make  util/myNewCustomFit.C


How to
---------
To run the code, after compiling (see "Setup"), use the command::

    ./myFit.exe  <action(s)>  [<config file>]  [<update>]  [<options>]

The configuration file (<config file>) is a text file containing all the information on the definition of samples and fit regions, including all the fit and draw options.
By default the file  config/myFit.config  is loaded.
See the section "Config File" for more details.
As examples give a look at  ``config/myFit.config``  or  ``config/ttH2015.config`` .
Most of the times the only file the user has to modify to get his fit is the configuration file.

The only mandatory argument, <action(s)>, tells to the TtHFitter which operation(s) to perform.
The possible operations are defined in the main file (e.g. util/myFit.C).
For instance, if you use the default file `util/myFit.C`, the available options are:

* h : read input histograms (valid only if the proper option is spefified in the config file)
* n : read input ntuples (valid only if the proper option is spefified in the config file)
* w : create the RooStats xmls and workspace
* f : fit the workspace
* l : calculate exclusion limit
* s : calculate significance
* d : draw pre-fit plots
* p : draw post-fit plots
* a : draw separation plots
* r : draw ranking plot (see later)
* b : re-run smoothing (in the future also rebinning)
* m : multi-fit (see later)

New optional argument: <options>.
It's a string (so make sure to use " or ' to enclose the string if you use more than one option) defining a list of options, in the form::

    "<option1>=<value1>,<value2>,...:<option2>=..."

See the section "Command line options" below.


Config File
---------

Here's a list of the inputs and options which can be specifed in the config file:

- The structure of the file should be the following::

     <ObjectType>: <ObjectName>
       <ObjectProperty>: <Value>
       <ObjectProperty>: <Value>
       ...

     <ObjectType>: <ObjectName>
       <ObjectProperty>: <Value>
       <ObjectProperty>: <Value>
       ...

     ...

NB: note the *blank* line between the objects!!)

- The file should contain
   * exactly one object of type "Fit"
   * at least one object of type "Sample"
   * at least one object of type "Region"
   * any number of objects of type "Systematic" (even 0 should be ok)

Note that, each object should have unique <ObjectName>.


- Then, for each object type, here's a PARTIAL list of properties to be specified:

  * Job:
     * Label: it's the label which will be showed on the plots
     * POI: the name of the parameter of interest; this should correspond to a NormFactor defined below
     * ReadFrom         : can be HIST or NTUP; default is HIST
     * HistoPath        : valid only for option HIST above is selected; it's the path where the input root files containing the histograms are stored
     * NtuplePath(s)    : valid only for option NTUP; it's the path(s) where the input root files containing the ntuples are stored
     * MCweight         : only for option NTUP; string defining the weight (for MC samples only)
     * Selection        : only for option NTUP; string defining the selection
     * NtupleName       : default name of the tree
     * Lumi             : value to scale all the "NormalizedByTheory" samples
     * LumiScale        : additional value to scale 'after' histogram creation (for fast scaling) IMPORTANT: use it only if you know what you are doing!!
     * SystPruningShape : Lower threshold to remove a shape systematic from the fit/limit (suppression is done per sample and per region) (Ex: 0.02 for 2%)
     * SystPruningNorm  : Lower threshold to remove a normalisation systematic from the fit/limit (suppression is done per sample and per region) (Ex: 0.02 for 2%)
     * SystLarge        : all systematics above this threshold will be flagged in the pruning plot) (e.g. 0.4 will flag systematics that are larger than 40%)
     * IntCodeOverall   : interpolation code used for the normalization component of systematics (should match the one used in RooStats)
     * IntCodeShape     : interpolation code used for the shape component of systematics (should match the one used in RooStats)
     * MCstatThreshold  : if set it will add the MC stat uncertainty to the fit (and to the plots); a NP will be added for each bin with an MC stat uncertainty > this threshold (relative)
     * DebugLevel       : 0 or 1
     * PlotOptions      : a set of options for plotting:
        * YIELDS : if set, the legend will be one-column and will include the yileds; otherwise two-columns and no yields
        * NORMSIG : add normlised signal to plots
        * NOSIG: don't show signal in stack
        * OVERSIG: overlay signal (not normalised)
        * CHI2: the chi2/ndf and chi2 prob will be printed on each plot, provided that the option GetChi2 is set
     * TableOptions      : a set of options for tables (yield tables only for the moment):
        * STANDALONE : default! If not set, no "\begin{document}"
        * FOOTNOTESIZE : -> \footnotesize
        * LANDSCAPE : -> \begin{landscape}
     * SystControlPlots : if set to true, plots will be dumped showing the shape effect of a given systematic (before and after smoothing/symmetrisation)
     * SystDataPlots    : if set to true, plots will be dumped showing the shape effect of a given systematic (before and after smoothing/symmetrisation) on top of the nominal sum of samples.Data are then plotted in the ratio. If the option is set to "fillUpFrame" data will also be plotted in the upper frame.
     * CorrelationThreshold : Threshold used to draw the correaltion matrix (only systematics with at least one correlation larger than than draw) (0.05:5%)
     * SignalRegionsPlot: list of regions to put in SignalRegionsPlot; use "EMPTY" to put an empty entry, "ENDL" to specify end of line
     * HistoChecks      : NOCRASH: means that if an error is found in the input histograms, the code continues (with only warnings) -- default leads to a crash in case of problem
     * LumiLabel        : label for luminosity to be put on plots
     * CmeLabel         : label for center-of-mass energy to be put on plots
     * SplitHistoFiles  : set this to true to have histogram files split by region (useful with many regions and/or run in parallel)
     * BlindingThreshold: bins with S/B > this number will be blinded
     * KeepPrefitBlindedBins : if set to TRUE, and if pre-fit an post-fit plots are produced together ("dp" opiton) pre-fit blinding is kept in post-fit plots 
     * RankingMaxNP     : max number of NP to show in ranking plot
     * RankingPlot      : NP categories in gammas or systs, if set to Systs(Gammas) then plot only systs(Gammas) in ranking, default produce plot for systs+gammas, can also set to all to have the 3 plots.
     * ImageFormat      : png, pdf or eps
     * StatOnly         : the code ignores systematics and MC stat uncertainties from all computations (limits, significances, fit, ...); need to re-reate ws in case of limit and sigificance
     * SystErrorBars    : add stat error bars to syst variations in syst plots
     * SummaryPlotRegions : list of regions to be shown in summary plot (useful to specify a custom order)
     * FixNPforStatOnly : if set to TRUE, when running stat-only (with either of the two options) also the norm facotrs other than the POI are kept fixed
     * InputFolder      : specify it to read fit input histograms from a different direcotry than <jobName>/Histograms/
     * InputName        : specify it to read fit input histograms from files with different name than <jobName>_blabla.root
     * OutputDir        : specify it to write everything in a different direcotry than <jobName>
     * WorkspaceFileName : if specified, an external ws can be used as input for fitting (not 100% supported)
     * KeepPruning      : if set to TRUE, the first time the ws is created (option w) a Pruning.root file is created under <jobName>/ and used for future operations to skip pruned systematics (makes operations much faster in case many syst are pruned)
     * AtlasLabel       : to specify Internal, Preliminary, etc...
     * CleanTables      : if set to TRUE, a cleaned version of the tex tables is created (basically removing the "#") - to be expanded
     * SystCategoryTables : if set to TRUE, additional syst tables with systematics grouped by category are created
     * SummaryPlotYmax  : if set, it will force the summary plot to use this value as max y-maxis value
     * SummaryPlotYmin  : if set, it will force the summary plot to use this value as min y-maxis value
     * RatioYmax        : if set, it will specify a the max of the range of the ratio plots
     * RatioYmin        : if set, it will specify a the min of the range of the ratio plots 
     * RatioYmaxPostFit : if set, it will specify a the max of the range of the ratio plots, for post-fit only
     * RatioYminPostFit : if set, it will specify a the min of the range of the ratio plots, for post-fit only 
     * CustomAsimov     : if set, the workspace will be created with an AsimovData built according to Sample->AsimovReplacementFor option (see below) instead of data
     * RandomPOISeed    : if set to a >= 0 number, the signal sample(s) to which the POI is assigned get scaled by a random number generated statrging from this seed, just before the ws creation; if the same seed is used in the cofig, post-fit plots will show consistent results (i.e. before post-fit drawing the POI is scaled by the same number)
     * GetChi2          : if set to TRUE (or STAT+SYST), for pre- and post-fit plots the extended chi2 test is done, and results are printed on the screen for each plot when running d and/or p; can be set to STAT (or STAT-ONLY) for stat-only chi2
     
  * Fit:
     * FitType          : can be SPLUSB (default) or BONLY to fit under the s+b or the b-only hypothesis
     * FitRegion        : can be CRSR (default) or CRONLY to fit considering both signal and control regions in the fit, or only control regions. You can also specify a coma-separated list of regions to use in the fit
     * FitBlind         : specify is real data or Asimov data should be used in the fit (TRUE or FALSE). By default, fit are NOT blind.
     * POIAsimov        : value of the parameter of interest in the AsimovDataset used in the fit
     * NPValues         : values of the nuisance parameters used to build the fit. Coma-separated list of NP:value (e.g. alpha_ttbarbb_XS:1,alpha_ttbarbcc_XS:1.5)
     * FixNPs           : values of the nuisance parameters used to be fixed in the fit. Coma-separated list of NP:value (e.g. alpha_ttbarbb_XS:1,alpha_ttbarbcc_XS:1.5)
     * doLHscan         : coma separated list of names of the POI or NP from which you want to produce the likelihood scan, if first element of the list is "all" then all systematics are profiled
     * UseMinos         : coma separated list of names of the POI and/or NP for which you want to calculate the MINOS errors, if first element of the list is "all" then the MINOS errors is calculated for all systematics and POIs
     * SetRandomInitialNPval : useful to set this to >0 (e.g. 0.1) to help convergence of Asimov fits
     * NumCPU           : specify the number of CPU to use for the minimization (default = 1)
     * StatOnlyFit      : if specified, the fit will keep fixed all the NP to the latest fit result, and the fit results will be saved with the _statOnly suffix (also possible to use it from command line)
     * GetGoodnessOfFit : set to TRUE to get it (based on chi2 probability from comparison of negative-log-likelihoods)
     
  * Limit:
     * LimitType        : can be ASYMPTOTIC or TOYS (the latter is not yet supported)
     * LimitBlind       : can be TRUE or FALSE (TRUE means that ALL regions are blinded)
     * POIAsimov        : value of the POI to inject in the Asimov dataset in LimitBlind is set to TRUE
     * SignalInjection  : if set to TRUE, expected signal with signal injection is evaluated

  * Options:
     * additional options, accepting only float as arguments - useful for adding your functionalities & flags in a quick way, since they need minimal changes in the code) ...

  * Region:
     * VariableTitle    : it's the label which will be displayed on the x-axis in the plots
     * Label            : it's the label which will be showed on the plots and specifies which region is shown
     * TexLabel         : label for tex files
     * ShortLabel       : same as above, but a shorter version for plots with smaller available place
     * LumiLabel        : label for luminosity to be put on plots
     * CmeLabel         : label for center-of-mass energy to be put on plots
     * LogScale         : set it to true to have log-scale when plotting this region
     * HistoFile        : only for option HIST, the file name to be used
     * HistoName        : only for option HIST, the histogram name to be used
     * HistoPathSuff(s) : only for option HIST, the path suffix (or suffixes, comma-separated) where to find the histogram files for this region
     * Variable         : only for option NTUP, the variable (or expression) inside the ntuple to plot can define a variable as X|Y to do the correlation plot between X and Y
     * VariableForSample: only for option NTUP, allows to set exceptions for Variable. This is a very useful feature when using TRF only in some samples. Comma-separated list of sample:variable (e.g. wjets:met_met/1e3,zjets:Mbbb/1e).
     * Selection        : only for option NTUP, the selection done on the ntuple for this region
     * NtupleName       : only for option NTUP, the name of the tree for this region
     * NtuplePathSuff(s): only for option NTUP, the path sufix (or suffixes, comma-separated) where to find the ntuple files for this region
     * MCweight         : only for option NTUP, the additional weight sed in this region (for MC samples only)
     * Rebin            : if specified, the histograms will be rebinned merging N bins together, where N is the argument (int)
     * Binning          : if specified, the histograms will be rebinned according to the new binning specifed, in the form like (0,10,20,50,100). If option AutoBin is set, use algorithms/functions ro define the binning. Example - Binning: "AutoBin","TransfoD",5.,6. (TransfoF also available, 5. and 6. are parameters of the transformation). If used in background region and zSig!=0 (first parameter, =0 gives flat background) then need a coma separated list of backgrounds to use instead of signal to compute the binning.
     * BinWidth         : if specified, two things are done: this number is used to decorate the y axis label and the bin content is scaled for bins with a bin width different from this number
     * Type             : can be SIGNAL, CONTROL or VALIDATION; used depending on Fit->FitType; if VALIDATION is set, the region is never fitted; default is SIGNAL
     * DataType         : ASIMOV or DATA. Is asimov is set, the limits and significances are computed without taking into acount the data in these region, but a projection of the fit performed in the regions with DATA
     * Ymax             : if set, it will force the plot to use this value as max y-maxis value
     * Ymin             : if set, it will force the plot to use this value as min y-maxis value
     * RatioYmax        : if set, it will specify a the max of the range of the ratio plot for this region only
     * RatioYmin        : if set, it will specify a the min of the range of the ratio plot for this region only 
     * RatioYmaxPostFit : if set, it will specify a the max of the range of the ratio plot for this region only, for post-fit only
     * RatioYminPostFit : if set, it will specify a the min of the range of the ratio plot for this region only, for post-fit only 
     * DropBins         : allows to specify a comma-separated list of bins to set to 0 (both for data and prediction), starting from 0 for the index

  * Sample:
     * Type             : can be SIGNAL, BACKGROUND, DATA or GHOST; default is BACKGROUND; GHOST means: no syst, not drawn, not propagated to workspace
     * Title            : title shown on the legends
     * TexTitle         : title shown on tex tables
     * Group            : if specified, sample will be grouped with other samples with same group and this label will be used in plots
     * HistoFile        : valid only for option HIST; which root file to read (excluding the suffix ".root"); this will be combined with Fit->HistoPath to build the full path
     * HistoName        : valid only for option HIST; name of histogram to read
     * NtuplePath       : valid only for option HIST; it's the path where the input root files containing the histograms are stored
     * NtupleFile(s)    : valid only for option NTUP; it's the file name(s) where the input ntuples are stored
     * NtupleName       : valid only for option NTUP; name of tree to read
     * NtuplePath(s)    : valid only for option NTUP; it's the path(s) where the input root files containing the ntuples are stored
     * FillColor        : histogram fill color (not valid for data)
     * LineColor        : histogram line color
     * NormFactor       : NormalisationFactor (free parameter in the fit); in the format <name>,nominal,min,max
     * NormalizedByTheory: set it to false for data-driven backgrounds (MCweight, Lumi and LumiScale from Job and Region will be ignored)
     * MCweight         : only for option NTUP, the additional weight sed in this sample (for all types of samples!! Not only MC)
     * Selection        : valid only for option NTUP; additional selection for this region
     * Regions          : set this to have the sample only in some regions
     * Exclude          : set this to exclude the sample in some regions
     * LumiScale(s)     : set this to scale the sample by a number; if more numbers are set, use a different one for each file / name / path...
     * IgnoreSelection  : if set, selection from Job and Region will be ignored
     * UseMCstat        : if set to FALSE, makes the fitter ignore the stat uncertainty for this sample
     * MultiplyBy       : if specified, each sample hist is multiplied bin-by-bin by another sample hist, in each of the regions
     * DivideBy         : if specified, each sample hist is divided bin-by-bin by another sample hist, in each of the regions
     * Smooth           : if set to TRUE, the nominal histograms are smoothed (based on TH1::Smooth but taking into account the original stat uncertainty) 
     * AsimovReplacementFor: only for GHOST samples; if set, the creation of cutsom Asimov data-set(s) is triggered; use as 'AsimovReplacementFor: "dataset","sample"', where "dataset" is the name of a custom Asimov dataset one wants to create (the same name will have to be set under Job->CustomAsimov in order to use it) and "sample" is the sample this GHOST sample will superseed

  * NormFactor:
     * Samples          : comma-separated list of samples on which to apply the norm factor
     * Regions          : comma-separated list of regions where to apply the norm factor
     * Exclude          : comma-separated list of samples/regions to exclude
     * Title            : title of the norm factor
     * Nominal          : nominal value
     * Min              : min value
     * Max              : max value
     * Constant         : set to TRUE to have a fixed norm factor

  * ShapeFactor:
     * Samples          : comma-separated list of samples on which to apply the shape factor
     * Regions          : comma-separated list of regions where to apply the shape factor
     * Title            : title of the shape factor

  * Systematic:
     * Samples             : comma-separated list of samples on which to apply the systematic
     * Regions             : comma-separated list of regions where to apply the systematic
     * Exclude             : comma-separated list of samples/regions to exclude
     * ExcludeRegionSample : comma-separated list of region:sample to exclude
     * Type                : can be HISTO, OVERALL or SHAPE (this refers to the HistFactory Shape Systematic, i.e. uncorrelated bin-by-bin)
     * Title               : title of the systematic (will be shown in plots)
     * Category            : major category to which the systematic belongs (instrumental, theory, ttbar, ...): used to split pulls plot for same category
     * HistoPathUp         : only for option HIST, for HISTO or SHAPE systematic: histogram file path for systematic up variation
     * HistoPathDown       : only for option HIST, for HISTO or SHAPE systematic: histogram file path for systematic down variation
     * HistoPathSufUp      : only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic up variation
     * HistoPathSufDown    : only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic down variation
     * HistoFileUp         : only for option HIST, for HISTO or SHAPE systematic: histogram file name for systematic up variation
     * HistoFileDown       : only for option HIST, for HISTO or SHAPE systematic: histogram file name for systematic down variation
     * HistoFileSufUp      : only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic up variation
     * HistoFileSufDown    : only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic down variation
     * HistoNameUp         : only for option HIST, for HISTO or SHAPE systematic: histogram name for systematic up variation
     * HistoNameDown       : only for option HIST, for HISTO or SHAPE systematic: histogram name for systematic down variation
     * HistoNameSufUp      : only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram names for systematic up variation
     * HistoNameSufDown    : only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram names for systematic down variation
     * NtuplePath(s)Up     : only for option NTUP, for HISTO or SHAPE systematic: ntuple file path(s) for systematic up variation
     * NtuplePath(s)Down   : only for option NTUP, for HISTO or SHAPE systematic: ntuple file path(s) for systematic down variation
     * NtuplePathSufUp     : only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file paths for systematic up variation
     * NtuplePathSufDown   : only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file paths for systematic down variation
     * NtupleFile(s)Up     : only for option NTUP, for HISTO or SHAPE systematic: ntuple file name(s) for systematic up variation
     * NtupleFile(s)Down   : only for option NTUP, for HISTO or SHAPE systematic: ntuple file name(s) for systematic down variation
     * NtupleFileSufUp     : only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file names for systematic up variation
     * NtupleFileSufDown   : only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file names for systematic down variation
     * NtupleName(s)Up     : only for option NTUP, for HISTO or SHAPE systematic: ntuple name(s) for systematic up variation
     * NtupleName(s)Down   : only for option NTUP, for HISTO or SHAPE systematic: ntuple name(s) for systematic down variation
     * NtupleNameSufUp     : only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple names for systematic up variation
     * NtupleNameSufDown   : only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple names for systematic down variation
     * WeightUp            : only for option NTUP, for HISTO or SHAPE systematic: weight for systematic up variation
     * WeightDown          : only for option NTUP, for HISTO or SHAPE systematic: weight for systematic down variation
     * WeightSufUp         : only for option NTUP, for HISTO or SHAPE systematic: additional weight for systematic up variation
     * WeightSufDown       : only for option NTUP, for HISTO or SHAPE systematic: additional weight for systematic down variation
     * IgnoreWeight        : only for option NTUP: if set, the corresponding weight (present in Job, Sample or Region) will be ignored for this systematic
     * Symmetrisation      : can be ONESIDED or TWOSIDED (...); for no symmetrisation, skip the line
     * Smoothing           : smoothing code to apply; use 40 for default smoothing; for no smoothing, skip the line
     * OverallUp           : for OVERALL systematic: the relative "up" shift (0.1 means +10%)
     * OverallDown         : for OVERALL systematic: the relative "down" shift (-0.1 means -10%)
     * ScaleUp             : for OVERALL, HISTO or SHAPE systematic: scale difference between "up" and nominal by a factor
     * ScaleDown           : for OVERALL, HISTO or SHAPE systematic: scale difference between "down" and nominal by a factor
     * ReferenceSample     : if this is specified, the syst variation is evaluated w.r.t. this reference sample (often a GHOST sample) instead of the nominal, and then the relative difference is propagated to nominal; NOTE: also the overall relative difference is propagated
     * KeepNormForSamples  : list of samples (or sum of samples, in the form smp1+smp2), comma separated, for which the systematic gets shape only in each region
     * PreSmoothing        : if set to TRUE, a TH1::Smooth-based smoothing is applied, prior to the usual smoothing (if set)
     * SubtractRefSampleVar: if set to TRUE, the relative variation of the ReferenceSample will be linearily subtracted from the relative variation of each affected sample, for the same systematic - this is relevant e.g. for Full JER SmearingModel, where data would be the reference sample
     

Command line options
---------

Currently the supported options are:

* Regions:     to limit the regions to use to the list specified
* Samples:     to limit the samples to use to the list specified
* Systematics: to limit the systematics to use to the list specified
* Signal:      in case more than one SIGNAL sample is specified in your config file, you can specify which one you want to run on (for plots, workspace creation and fits/limits/significance)
* Exclude:     to exclude certain Regions / Samples / Systematics
* Suffix:      used for: plots, workspace, fit resutls, etc
* SaveSuffix:  used for: saving histograms with a suffux (to be merged / renamed later, see last section on hupdate)
* Update:      if TRUE, the output .root file is updated, otherwise is overwrote
* StatOnlyFit: if TRUE, the same as Fit->StatOnlyFit

Note: the wild-card * is supported, but only as last character.
Example::

     ./myFit.exe  n  config/ttH2015.config 'Regions=HThad_ge6jge4b:Exclude=BTag_*'


Ranking Plot
---------

- The ranking plot can be created in one go, with just the command line argument "r" (after having run the nominal fit fit "f").
- Since this can take too much time (and memory), for complicated fits it's better to run it in several steps:
   by specifying the command-line option "Ranking=<name/index>", one can produce the txt input for the ranking only for a specific line of the ranking, i.e. for a single NP (speficied either through its name or index). Once all the needed txt files are created (e.g. in parallel through batch jobs) with the option "Ranking=plot" they are merged to create the final plot.

- Examples::

     # this runs the ranking in one go
     ./myFit.exe  r  <config>
     #these commands will first create the inputs for the ranking one by one and then merge them in the plot
     ./myFit.exe  r  <config> Ranking=Lumi
     ./myFit.exe  r  <config> Ranking=JES1
     ./myFit.exe  r  <config> Ranking=ttXsec
     ./myFit.exe  r  <config> Ranking=plot


Multi-Fit
---------

The Multi-Fit functionality can be used to compare fit results or even to combine fit inputs from different configuration files / Jobs.

- To use it you need a dedicated config file, with a similar structure as the usual ones. Example::

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

- This config file can be run with the command line::

    ./myFit  m  config/myTopWS_multifit.config

  this will compare the fit resutls in terms of fitted NP, fitted POI and limits from the two config files specified. Notice that the fit and limits results have to be already available (they are not produced on the flight).

- To make a real combination, one needs to use the usual command options "w", "f" and "l" together with the flag "Combine: TRUE" in the config above. Example::

    ./myFit  mwf  config/myTopWS_multifit.config

  this will create a combined ws starting from the individual ws for the different regions in the two config files, and fit it.


Input File Merging with hupdate
---------

- A macro hupdate is included in order to mimic hadd functionality but without adding histograms if they have the same name
- This is useful for running different systematics in different steps (like different batch jobs) and then merging results afterward
- One needs to compile the hupdate macro before being able to use it::

    make hupdate

- Example usage, combined with the usage of SaveSuffix::

    make hupdate
    ./myFit.exe n config/ttH2015.config Systematics=BTag_B_NP1:SaveSuffix=_BTag_B_NP1
    ./myFit.exe n config/ttH2015.config Exclude=BTag_B_NP1:SaveSuffix=_rest
    ./hupdate ttH2015/Histograms/ttH2015_HThad_4j2b_histos.root ttH2015/Histograms/ttH2015_HThad_4j2b_histos_rest.root ttH2015/Histograms/ttH2015_HThad_4j2b_histos_BTag_B_NP1.root
    ./hupdate ttH2015/Histograms/ttH2015_HThad_5j3b_histos_NEW.root ttH2015/Histograms/ttH2015_HThad_5j3b_histos.root ttH2015/Histograms/ttH2015_HThad_5j3b_histos_BTag_B_NP1.root
    ./hupdate ttH2015/Histograms/ttH2015_HThad_ge6jge4b_histos_NEW.root ttH2015/Histograms/ttH2015_HThad_ge6jge4b_histos.root ttH2015/Histograms/ttH2015_HThad_ge6jge4b_histos_BTag_B_NP1.root
    ./myFit.exe dwf config/ttH2015.config


Output Directories Structure
---------

* For each TtHFit objetc, a diretory is created, with the same name as the Fit Name
* Inside this direcotry, at every step, some outputs are created, following the structure described above

  * Plots/              : contains the data/MC plots, pre- and post-fit, for all the Signal, Control and Validation regions, including the summary plots
  * Tables/             : contains the tables in txt and tex format
  * RooStats/           : contains the workspace(s) and the xmls
  * Fits/               : contains the output from fits
  * Limits/             : contains the outputs from the limit-setting code
  * Significance/       : contains the outputs from the significance code
  * Systematics/        : contains the plots for the syst variations
  * Histograms/         : contains the root file(s) with all the inputs
  * LHoodPlots/         : contains the likelihood scan with respect to the specified parameter

ShapeFactor example
-------------------

* The following scripts create example histograms in :code:`exampleDataDriven` directory and execute  :code:`myFit.exe` using :code:`config/dataDriven.config`
* The example contains a control region and signal region with two bins. The shape of one of the background samples is estimated using the ShapeFactor::

    python makeDataDriven.py
    python runDataDrivenExample.py


The results are in :code:`JobDataDriven`


TtHFitter package authors
-----------------
Contacts:

* Michele Pinamonti <michele.pinamonti@gmail.com>
* Loic Valery <loic.valery@cern.ch>
