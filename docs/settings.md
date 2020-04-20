# Settings for TRExFitter blocks

## Standard fit
The following settings are for normal fits, performed without the action `m`.

### `Job` block settings
| **Option** | **Function** |
| ---------- | ------------ |
| Label                        | the label which will be shown on plots, if 'none' is set, no label will be shown |
| POI                          | the name of the parameter of interest; this should correspond to a NormFactor defined below |
| ReadFrom                     | can be HIST or NTUP; default is HIST |
| HistoPath(s)                 | valid only for option HIST above is selected; it's the path(s) where the input root files containing the histograms are stored |
| HistoFile(s)                 | valid only for option HIST; it's the file name(s) where the input root files containing the histograms are stored |
| HistoName(s)                 | valid only for option HIST; it's the histogram name(s) to read from the file(s) |
| HistoNameNominal             | valid only for option HIST; name of the nominal histogram, in case the systematic histogram names cannot be build by suffixing but instead replace the nominal name. Behaves like HistoNameSuff for Systematics if Nominal was a Systematic |
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
| PruningShapeOption           | Can be `MAXBIN` (Default) or `KSTEST`. `MAXBIN` checks if the systematic variation has at least one bin which is different compared to the nominal (variation is >= threshold). `KSTEST` runs the Kolmogorov test (using the `X` option, a.k.a running pseudoexperiments). When using `KSTEST`, if the probability is <= (1-threshold) then the systematic is marked as having shape component. |
| DoSystNormalizationPlots     | Set to `FALSE` to disable the normalization summary plot that is produced during the `w` step |
| IntCodeOverall               | interpolation code used for the normalization component of systematics (should match the one used in RooStats) |
| IntCodeShape                 | interpolation code used for the shape component of systematics (should match the one used in RooStats) |
| MCstatThreshold              | by default, the MC stat uncertainty is included in the fit (and to the plots); a NP will be added for each bin with an MC stat uncertainty > this threshold (relative) if the option is set to a float (default: no threshold); can also set to NONE in order to disable MC stat uncertainty completely |
| MCstatConstraint             | constraint used for MC stat uncertainties, can be set to 'GAUSSIAN' or 'POISSON' (default) |
| DebugLevel                   | 0 = prints only Warning and Errors, 1 = additionally prints Info messages, 2 = additionally prints Debug messages, >2 additionally prints Verbose messages. For option <2 RooFit/Roostats messages will be heavily suppressed |
| Logo                         | is set to TRUE will print the TRExFitter logo |
| PlotOptions                  | a set of options for plotting:<br>&nbsp; &nbsp; **YIELDS**: if set, the legend will be one-column and will include the yields; otherwise two-columns and no yields<br>&nbsp; &nbsp; **NORMSIG**: add normlised signal to plots<br>&nbsp; &nbsp; **NOSIG**: don't show signal in stack<br>&nbsp; &nbsp; **OVERSIG**: overlay signal (not normalised)<br>&nbsp; &nbsp; **CHI2**: the chi2/ndf and chi2 prob will be printed on each plot, provided that the option GetChi2 is set<br>&nbsp; &nbsp; **PREFITONPOSTFIT**: draw a dashed line on the postfit plot that indicates the sum of prefit background<br>&nbsp; &nbsp; **NOXERR**: removes the horizontal error bars on the data and the ratio plots |
| PlotOptionsSummary           | the same as PlotOptions but for the summary plot (if nothing is specified, PlotOptions is used) |
| RatioYtitle                  | Label to be used on the Y-axis of the ratio plot, default is "Data/pred." |
| TableOptions                 | a set of options for tables:<br>&nbsp; &nbsp; **STANDALONE**: default! If not set, no "\begin{document}"<br>&nbsp; &nbsp; **FOOTNOTESIZE**: -> \footnotesize <br>&nbsp; &nbsp; **LANDSCAPE**: -> \begin{landscape} |
| SystControlPlots             | if set to TRUE, plots showing the shape effect of a given systematic (before and after smoothing/symmetrisation) will be produced (default: TRUE) |
| SystErrorBars                | TRUE by default, add stat error bars to syst variations in syst plots, set to FALSE to disable |
| SystDataPlots                | if set to TRUE, plots showing the shape effect of a given systematic (before and after smoothing/symmetrisation) on top of the nominal sum of samples will be produced. Data are then plotted in the ratio. If the option is set to "fillUpFrame", data will also be plotted in the upper frame. |
| CorrelationThreshold         | Threshold used to draw the correlation matrix (only systematics with at least one correlation larger than than draw) (0.05:5%) |
| SignalRegionsPlot            | list of regions to put in SignalRegionsPlot and PieChartPlots; use "EMPTY" to put an empty entry, "ENDL" to specify end of line. This specifies the order of regions plotted in signal region S/B plots and pie chart plots, as well as number of regions per row. |
| HistoChecks                  | NOCRASH: means that if an error is found in the input histograms, the code continues (with only warnings) -- default leads to a crash in case of problem, if set to NOCRASH, also prints warning instead of error (and crash) when input files are not found for the histogram building step |
| LumiLabel                    | label for luminosity to be put on plots |
| CmeLabel                     | label for center-of-mass energy to be put on plots |
| SplitHistoFiles              | set this to TRUE to have histogram files split by region (useful with many regions and/or run in parallel) |
| BlindingThreshold            | blind bins when S/B is greater than this threshold, use the `BlindingType` option for other definitions than just S/B |
| BlindingType                 | how to calculate the quantity to determine blinding, options are SOVERB (for S/B), SOVERSPLUSB(for S/(S+B)), SOVERSQRTB (for S/sqrt(B)) and SOVERSQRTSPLUSB (for S/sqrt(S+B)), default is SOVERB |
| KeepPrefitBlindedBins        | if set to TRUE, and if pre-fit an post-fit plots are produced together ("dp" option) pre-fit blinding is kept in post-fit plots |
| RankingMaxNP                 | max number of NP to show in ranking plot |
| RankingPlot                  | NP categories in gammas or systs, if set to Systs(Gammas) then plot only systs(Gammas) in ranking, default produce plot for systs+gammas, can also set to all to have the 3 plots. |
| ImageFormat                  | png, pdf or eps |
| StatOnly                     | the code ignores systematics and MC stat uncertainties from all computations (limits, significances, fit, ...); need to re-create ws in case of limit and significance |
| SummaryPlotRegions           | list of regions to be shown in summary plot (useful to specify a custom order) |
| FixNPforStatOnly             | if set to TRUE, when running stat-only (with either of the two options) also the norm factors other than the POI are kept fixed |
| InputFolder                  | specify it to read fit input histograms from a different directory than `<jobName>/Histograms/` |
| InputName                    | specify it to read fit input histograms from files with different name than `<jobName>_blabla.root` |
| OutputDir                    | specify it to write everything in a different directory, the full directory structure will be: "OutputDirName/" (if set) + jobname (can be set via a command line option 'Job') |
| WorkspaceFileName            | if specified, an external ws can be used as input for fitting (not 100% supported) |
| KeepPruning                  | if set to TRUE, the first time the ws is created (option w) a Pruning.root file is created under `<jobName>/` and used for future operations to skip pruned systematics (makes operations much faster in case many syst are pruned) |
| AtlasLabel                   | to specify Internal, Preliminary, etc... If set to `none` the whole label will be removed |
| CleanTables                  | if set to TRUE, a cleaned version of the tex tables is created (basically removing the "#") - to be expanded |
| SystCategoryTables           | if set to TRUE, additional syst tables with systematics grouped by category are created |
| SummaryPlotYmax              | if set, it will force the summary plot to use this value as max y-maxis value |
| SummaryPlotYmin              | if set, it will force the summary plot to use this value as min y-maxis value |
| SummaryLogY                  | Can be set to `TRUE` (default) or `FALSE`, controls the use of logarithmic scale on Y axis for summary plots |
| RatioYmax                    | if set, it will specify the max of the range of the ratio plots |
| RatioYmin                    | if set, it will specify the min of the range of the ratio plots |
| RatioYmaxPostFit             | if set, it will specify the max of the range of the ratio plots, for post-fit only |
| RatioYminPostFit             | if set, it will specify the min of the range of the ratio plots, for post-fit only |
| CustomAsimov                 | if set, the workspace will be created with an AsimovData built according to Sample->`AsimovReplacementFor` option (see below) instead of data |
| GetChi2                      | if set to TRUE (or STAT+SYST), for pre- and post-fit plots the extended chi2 test is done, and results are printed on the screen for each plot when running d and/or p; can be set to STAT (or STAT-ONLY) for stat-only chi2 |
| SmoothingOption              | Choose which smoothing option to use, allowed parameters are: MAXVARIATION (default), TTBARRESONANCE (see also [FAQ section](#faq)), COMMONTOOLSMOOTHMONOTONIC, COMMONTOOLSMOOTHPARABOLIC, KERNELRATIOUNIFORM, KERNELDELTAGAUSS or KERNELRATIOGAUSS. |
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
| Bootstrap                    | (only works with NTUP inputs) if set, the bootstrap method wil be used; the argument should be a string like `bsWeight(BootstrapIdx,eventNumber,mcChannelNumber)`, where `bsWeight` should be loaded with `CustomFunctions: "bsWeight.C"` and eventNumber and mcChannelNumber should be existing branches for all the MC ntuples; then, to produce the i-th bootstrap pseudo-experiment, or to run on it (e.g. to perform a fit) the command-line option `BootstrapIdx=<i>` should be given, with `<i>=0,1,2,3...`. Alternatively there is default function `util/BootstrapDefault.C` which would be called as `BootstrapDefault(eventNumber+mcChannelNumber+BootstrapIdx)`. |
| BootstrapSyst                | Specifies the systematics for which should be the bootstrap done. |
| BootstrapSample              | Specifies the sample for which should be the bootstrap done. Nominal and all correlated systematics are smeared, correlated systematics are specificied using IsCorrelated option. |
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
| SmoothMorphingTemplates      | if set to TRUE (default is FALSE), the templates used for morphing are forced to have linear dependence on the morphing parameter, bin-by-bin (plots are produced per bin, in the Morphing directory) |
| PropagateSystsForMorphing    | if set to TRUE (default is FALSE), the non-nominal templates inherit all the systematics from the nominal template; NB: the nominal template is determined by the nominal value of the norm factor in the config |
| SummaryPrefix                | adds a prefix to summary and merge plots |
| AllowWrongRegionSample       | Can be TRUE or FALSE (default). When set to TRUE code will print only warnings when chosen samples or regions for various options are not defined. When set to FALSE the code will print errors and stop when the samples/regions are not defined. |
| POIPrecision                 | Integer value N, N >=1 and N <=5. Will tell the code to use N decimal places for norm factor mean value and uncertainty. Default is 2 |
| RankingPOIName               | Custom name for the POI for ranking plots. Default is `#mu` |
| UseGammasForCorr             | If set to `TRUE` will add gammas into correlation matrix plot. Default is `FALSE` |
| UseATLASRounding             | If set to `TRUE` will use PGD/ATLAS rounding to yield tables (both .txt and .tex) |
| UseATLASRoundingTxt          | If set to `TRUE` will use PGD/ATLAS rounding to yield tables (only .txt) |
| UseATLASRoundingTex          | If set to `TRUE` will use PGD/ATLAS rounding to yield tables (only .tex) |
| PrePostFitCanvasSize         | Set width and height for canvas for pre/post-fit plots  |
| SummaryCanvasSize            | Set width and height for canvas for summary plots  |
| MergeCanvasSize              | Set width and height for canvas for merged plots  |
| PieChartCanvasSize           | Set width and height for canvas for pie chart plots  |
| NPRankingCanvasSize          | Set width and height for canvas for NP ranking plot  |
| PruningType                 | Can be set to `BACKGROUNDREFERENCE` or `COMBINEDREFERENCE` (default is `SEPARATESAMPLE`), and pruning (both shape and norm) will be done w.r.t. to total/total-background. |
| LabelX                       | Custom X position for ATLAS label and others on Data/MC plots. |
| LabelY                       | Custom Y position for ATLAS label and others on Data/MC plots. |
| LegendX1                     | Custom Legend X1 position for ATLAS label and others on Data/MC plots. |
| LegendX2                     | Custom Legend X2 position for ATLAS label and others on Data/MC plots. |
| LegendY                      | Custom Legend Y top position for ATLAS label and others on Data/MC plots. |
| LabelXSummary                | Same as LabelX but for Summary plot. |
| LabelYSummary                | Same as LabelY but for Summary plot. |
| LegendX1Summary              | Same as LegendX1 but for Summary plot. |
| LegendX2Summary              | Same as LegendX2 but for Summary plot. |
| LegendYSummary               | Same as LegendY but for Summary plot. |
| LabelXMerge                  | Same as LabelX but for Merged plot. |
| LabelYMerge                  | Same as LabelY but for Merged plot. |
| LegendX1Merge                | Same as LegendX1 but for Merged plot. |
| LegendX2Merge                | Same as LegendX2 but for Merged plot. |
| LegendYMerge                 | Same as LegendY but for Merged plot. |
| LegendNColumns               | Number of columns in Legend for Data/MC plots. |
| LegendNColumnsSummary        | Same as LegendNColumns but for Summary plot. |
| LegendNColumnsMerge          | Same as LegendNColumns but for Merged plot. |
| ExcludeFromMorphing          | The specified sample is left our from the morphing (useful to do closure tests for morphing). |
| ScaleSamplesToData           | The specified samples will be scaled to data (when doing the d step). |
| MaxNtupleEvents              | valid only for option NTUP; if set to N, only first N entries per ntuple read (useful for debugging) |
| ResponseMatrixName(s)        | Name(s) of the histogram for response matrix |
| ResponseMatrixFile(s)        | File path(s) of the histogram for response matrix |
| ResponseMatrixPaths(s)       | Folder path(s) of the histogram for response matrix |
| ResponseMatrixNameNominal    | Nominal histogram name |
| AcceptanceName(s)        | Name(s) of the histogram for acceptance |
| AcceptanceFile(s)        | File path(s) of the histogram for acceptance |
| AcceptancePaths(s)       | Folder path(s) of the histogram for acceptance |
| AcceptanceNameNominal    | Nominal histogram name |
| SelectionEffName(s)       | Name(s) of the histogram for selection efficiency |
| SelectionEffFile(s)       | File path(s) of the histogram for selection efficiency |
| SelectionEffPaths(s)      | Folder path(s) of the histogram for selection efficiency |
| SelectionEffNameNominal   | Nominal histogram name |
| MigrationName(s)        | Name(s) of the histogram for migration matrix |
| MigrationFile(s)        | File path(s) of the histogram for migration matrix |
| MigrationPaths(s)       | Folder path(s) of the histogram for migration matrix |
| MigrationNameNominal    | Nominal histogram name |
| ReorderNPs              | If set to TRUE, fit results will show NPs and norm factors ordered according to how they appear in the config file. |
| BlindSRs                | If set, all SRs are forced to have DataType = ASIMOV |


### `Fit` block settings
| **Option** | **Function** |
| ---------- | ------------ |
| FitType                      | can be SPLUSB (default) or BONLY to fit under the s+b or the b-only hypothesis. Use UNFOLDING for unfolding|
| FitRegion                    | can be CRSR (default) or CRONLY to fit considering both signal and control regions in the fit, or only control regions. You can also specify a comma-separated list of regions to use in the fit |
| FitBlind                     | specify is real data or Asimov data should be used in the fit (TRUE or FALSE). By default, fit are NOT blind. |
| POIAsimov                    | value of the parameter of interest in the AsimovDataset used in the fit |
| NPValues                     | values of the nuisance parameters used to build the Asimov. Coma-separated list of NP:value (e.g. alpha_ttbarbb_XS:1,alpha_ttbarbcc_XS:1.5) |
| FixNPs                       | values of the nuisance parameters used to be fixed in the fit. Coma-separated list of NP:value (e.g. alpha_ttbarbb_XS:1,alpha_ttbarbcc_XS:1.5), currently only implemented for the `f` step |
| doLHscan                     | comma separated list of names of the POI or NP from which you want to produce the likelihood scan, if first element of the list is "all" then all systematics are profiled |
| do2DLHscan                   | produces 2D likelihood scan between the chosen parameters. Syntax: "paramX1,paramY1:param X2,paramY2". Warning takes long time. You can reduce the number of steps via `LHscanSteps`. Alternatively you can split up the 2D scan in slices with `Parallel2Dscan` |
| LHscanMin                    | minimum value for the LH scan on x-axis (default is Norm min). This also effect the x-axis in a 2D scan |
| LHscanMax                    | maximum value for the LH scan on x-axis (default is Norm max). This also effect the x-axis in a 2D scan |
| LHscanSteps                  | number of steps on the LH scan (default is 30) . This also effect the x-axis in a 2D scan ||
| LHscanMinY                   | minimum value for the 2D-LH scan on y-axis (default it Norm min) |
| LHscanMaxY                   | maximum value for the 2D-LH scan on y-axis (default is Norm max) |
| LHscanStepsY                 | number of steps on the LH scan in y-direction (default is 30, but if not specified it uses LHscanSteps) |
| UseMinos                     | comma separated list of names of the POI and/or NP for which you want to calculate the MINOS errors, if first element of the list is "all" then the MINOS errors is calculated for all systematics and POIs |
| SetRandomInitialNPval        | useful to set this to >0 (e.g. 0.1) to help convergence of Asimov fits |
| SetRandomInitialNPvalSeed    | seed used to determine initial NP settings in minimization process if SetRandomInitialNPval option is enabled |
| NumCPU                       | specify the number of CPU to use for the minimization (default = 1) |
| StatOnlyFit                  | if specified, the fit will keep fixed all the NP to the latest fit result, and the fit results will be saved with the `_statOnly` suffix (also possible to use it from command line) |
| GetGoodnessOfFit             | set to TRUE to get it (based on chi2 probability from comparison of negative-log-likelihoods) |
| SaturatedModel               | set it to TRUE to be able to get the goodness-of-fit test using the saturated model; if set to TRUE when running `w`, the resulting workspace will contain the saturated-model norm-factors; if set to TRUE when running `f` and `GetGoodnessOfFit` is set to TRUE as well, the goodness of fit is evaluated using the saturated model |
| DoNonProfileFit              | if set to TRUE (default is FALSE), instead of the fit profiling the systematics, a set of stat-only fits will be performed, on an Asimov data-set created with one syst variation at a time |
| FitToys                      | if set to N > 0, N stat-ony toys are generated and fitted |
| ToysHistoMin                 | If FitToys is used, set minimum on the output toys histogram X axis |
| ToysHistoMax                 | If FitToys is used, set maximum on the output toys histogram X axis |
| ToysHistoNbins               | If FitToys is used, set number of bins for toys histogram output |
| ToysPseudodataNP             | Name of the NP to be varied as pseudodata. Need to contain "alpha_NP" for NP called "NP". |
| ToysPseudodataNPShift        | Value of the NP to be used for pseudodata creation with "fToysPseudodataNP". Default value is 1 (represents pre-fit shift). |
| TemplateInterpolationOption  | Option only for morphing, tells the code which interpolation between the templates is used. Three possible options are available: LINEAR(default)/SMOOTHLINEAR/SQUAREROOT. All of these options basically use linear interpolation but SMOOTHLINEAR approximates it by integral of hyperbolic tangent and SQUAREROOT approximates it by $\sqrt{x^2+\epsilon}$ to achieve smooth transitions (first derivative) between the templates |
| BlindedParameters            | A comma separated list of POI/NPs that will be written as a hexadecimal number so it is not easy to read to not accidentally unblind. When at least one parameter is set the console output of the minimization is removed.
| DoNonProfileFitSystThreshold | When performing a NonProfileFit, systematics are not added to total if smaller than this threshold |
| NPValuesFromFitResults       | If set to a valid path pointing to a fit-result text file, the NPValues for Asimov-data creation will be readed from it |
| SetGlobalObservables         | If set to TRUE (default is FALASE), and if NPValues or NPValuesFromFitResults are set, also the global observables are shifted in the Likelihood according to the parameter values |


### `Limit block` settings
| **Option** | **Function** |
| ---------- | ------------ |
| LimitType                    | can be ASYMPTOTIC or TOYS (the latter is not yet supported) |
| LimitBlind                   | can be TRUE or FALSE (TRUE means that ALL regions are blinded) |
| POIAsimov                    | value of the POI to inject in the Asimov dataset in LimitBlind is set to TRUE |
| SignalInjection              | if set to TRUE, expected signal with signal injection is evaluated |
| SignalInjectionValue         | Value for the injected signal |
| ParamName                    | Name for the parameter in the output ROOT file |
| ParamValue                   | Value of the parameter in the output file (e.g. 172.5 for top mass) |
| OutputPrefixName             | Prefix for the output ROOT file |
| ConfidenceLevel              | Confidence level for the CLs. Default is 0.95 |


### `Significance` block settings
| **Option** | **Function** |
| ---------- | ------------ |
| SignificanceBlind            | can be TRUE or FALSE (TRUE means that ALL regions are blinded) |
| POIAsimov                    | value of the POI to inject in the Asimov dataset in SignificanceBlind is set to TRUE |
| ParamName                    | Name for the parameter in the output ROOT file |
| ParamValue                   | Value of the parameter in the output file (e.g. 172.5 for top mass) |
| OutputPrefixName             | Prefix for the output ROOT file |


### `Options` block settings
additional options, accepting only float as arguments - useful for adding your functionalities & flags in a quick way, since they need minimal changes in the code ...


### `Region` block settings
| **Option** | **Function** |
| ---------- | ------------ |
| VariableTitle                | it's the label which will be displayed on the x-axis in the plots |
| Label                        | it's the label which will be showed on the plots and specifies which region is shown |
| ResponseMatrixFile(s)        | same as for Job but for Region |
| ResponseMatrixName(s)        | same as for Job but for Region |
| ResponseMatrixPath(s)        | same as for Job but for Region |
| AcceptanceFile(s)            | same as for Job but for Region |
| AcceptanceName(s)            | same as for Job but for Region |
| AcceptancePath(s)            | same as for Job but for Region |
| SelectionEffFile(s)          | same as for Job but for Region |
| SelectionEffName(s)          | same as for Job but for Region |
| SelectionEffPath(s)          | same as for Job but for Region |
| MigrationEffFile(s)          | same as for Job but for Region |
| MigrationEffName(s)          | same as for Job but for Region |
| MigrationEffPath(s)          | same as for Job but for Region |
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
| NtuplePathSuff(s)            | only for option NTUP, the path suffix (or suffixes, comma-separated) where to find the ntuple files for this region |
| MCweight                     | only for option NTUP, the additional weight used in this region (for MC samples only) |
| Rebin                        | if specified, the histograms will be rebinned merging N bins together, where N is the argument (int) |
| Binning                      | if specified, the histograms will be binned according to the new binning specified, in the form like `0,10,20,50,100`. The parameter `AutoBin` allows to use algorithms to define the binning. Example - `Binning: "AutoBin","TransfoD",5.,6.` (`TransfoF` also available, `5.` and `6.` are parameters of the transformation; `TransfoJ` takes three parameters). A popular setting is `TransfoD` with both parameters set to the same value; the amount of bins used is equal to the sum of both parameters. If used in background region and zSig!=0 (first parameter, =0 gives flat background) then need a comma separated list of backgrounds to use instead of signal to compute the binning. More information about the automatic binning algorithms is found in the FAQ (see README). |
| Rebinning                    | if specified, the histograms will be rebinned according to the new binning specified, in the form like `0,10,20,50,100`. Differently from the BInning option, this one performs the rebinning aftre the original histograms are created. This means that this option can changed (or removed) before running the b step. |
| BinWidth                     | if specified, two things are done: this number is used to decorate the y axis label and the bin content is scaled for bins with a bin width different from this number |
| BinLabels                    | if specified, bin labels are set according to provided comma separated list (list length must be equal to number of bins) |
| Type                         | can be SIGNAL, CONTROL or VALIDATION; used depending on Fit->FitType; if VALIDATION is set, the region is never fitted; default is SIGNAL |
| DataType                     | ASIMOV or DATA. Is Asimov is set, limits, significance and final fit are computed without taking into account the data in this region, but creating a modified Asimov data-set through a projection of the fit performed in the regions with DATA only |
| Ymax                         | if set, it will force the plot to use this value as max y-maxis value |
| Ymin                         | if set, it will force the plot to use this value as min y-maxis value |
| RatioYmax                    | if set, it will specify the max of the range of the ratio plot for this region only |
| RatioYmin                    | if set, it will specify the min of the range of the ratio plot for this region only |
| RatioYmaxPostFit             | if set, it will specify the max of the range of the ratio plot for this region only, for post-fit only |
| RatioYminPostFit             | if set, it will specify the min of the range of the ratio plot for this region only, for post-fit only |
| DropBins                     | Allows to specify a comma-separated list of bins where the yield will be set to 0 (both for data and prediction), starting from 1 for the index as is the ROOT convention for bin indices. This can be used to remove bins from the fit, for example to fit a background-only model without bins that have significant signal contamination. |
| AutomaticDropBins            | Set this to TRUE to automatically drop bins (= set yield to 0 for data and prediction) above the `BlindingThreshold`. Default setting is FALSE, and manual use of `DropBins` also sets this option to FALSE. When FALSE, the bins are still blinded in plots via `BlindingThreshold`, but they will enter the fit (unless manually removed with `DropBins`). |
| Group                        | if specified, regions of the same group appear together in several places, see RegionGroups option |
| YaxisTitle                   | title of y-axis used for plots of the region |
| YmaxScale                    | scales range of y-axis (default: 2.0, meaning the maximum axis value is twice the largest yield in any bin) |
| Ymax                         | maximum value on y-axis |
| SkipSmoothing                | if smoothing of nominal samples is used, this option can be used to disable smoothing per region (default: FALSE) |
| XaxisRange                   | Manually call 'SetRangeUser()' on X axis. Needs two parameters(floats): min,max |
| NumberOfRecoBins             | Number of reco bins in this region when Unfolding is used |
| IsBinOfRegion                | can be set in order to declare this region as corresponding to a certain bin of another region (to be set as VALIDATION region and to be declared previously in the config file); in this way, one can perform a two-dimensional fit, by keeping the systematic smoothing in both directions: for each systematic, the smoothing applied to the VALIDATION region indicated will be propagated to the overall part of the systematic variation in this region; syntax: "region-name":binNmbr (NB: bin numbering here starts at 1 due to technical reasons) |


### `Sample` block settings
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
| FillColorRGB                 | histogram fill color in RGB (not valid for data). This expects a triplet of RGB values between 0 and 255, e.g. `255,0,0`. If set, the FillColor option is ignored. |
| LineColor                    | histogram line color |
| LineColorRGB                 | histogram line color in RGB. This expects a triplet of RGB values between 0 and 255, e.g. `255,0,0`. If set, the LineColor option is ignored. |
| NormFactor                   | NormalisationFactor (free parameter in the fit); in the format \<name\>,nominal,min,max |
| ShapeFactor                  | ShapeFactor added |
| NormalizedByTheory           | set it to FALSE for data-driven backgrounds (MCweight, Lumi and LumiScale from Job and Region will be ignored) |
| MCweight                     | only for option NTUP, the additional weight used in this sample |
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
| MCstatScale                  | scales up/down the MC stat size; useful to project sensitivity to larger MC samples |
| SystFromSample               | set it to TRUE (default FALSE) to have this sample inheriting all the systematics from another sample |


### `NormFactor` block settings
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
| Expression                   | a way to correlate this norm factor with other norm factors (using AddPreprocessFunction); two arguments, in the form `<expression>:<dependencies>`, where `<dependencies>` should contain the names of the norm factors the expression depends on, their nominal values and existence ranges [example: `(1.+Pmag*cos(theta))/2.:Pmag[0.9,0,1],theta[0,0,3.14]`] |
| Tau                          | If set, a constraint term will be added to the likelihood for this NF, in the same way as for the Tikhonov regularization (used in unfolding); the constraint will be of the form exp((1/2) * tau^2 * (NF-µ)^2), where µ is the Nominal value for the NF, so the larger is tau, the stronger the constraint |


### `ShapeFactor` block settings
| **Option** | **Function** |
| ---------- | ------------ |
| Samples                      | comma-separated list of samples on which to apply the shape factor |
| Regions                      | comma-separated list of regions where to apply the shape factor |
| Title                        | title of the shape factor |


### `Systematic` block settings
| **Option** | **Function** |
| ---------- | ------------ |
| Samples                      | comma-separated list of samples on which to apply the systematic |
| Regions                      | comma-separated list of regions where to apply the systematic |
| Exclude                      | comma-separated list of samples/regions to exclude |
| ExcludeRegionSample          | comma-separated list of region:sample to exclude |
| Type                         | can be HISTO, OVERALL, SHAPE (this refers to the HistFactory Shape Systematic, i.e. uncorrelated bin-by-bin) or STAT (this refers to auto-creation of one systematic from stat uncertainty for each bin of corresponding region - DEPRECATED) |
| Title                        | title of the systematic (will be shown in plots) |
| StoredName                   | if specified, will be used to read and write histograms in the root files under Histograms/ instead of the syst name; useful to decorrelate without re-creating histograms |
| NuisanceParameter            | if specified, this will be given to RooStats instead of the syst name; useful (and recommended) way to correlate systematics |
| IsFreeParameter              | if set to TRUE, the constraint will be a flat one instead of Gaussian (use with caution) |
| IsCorrelated                 | if set to TRUE, the systematic is considered correlated to nominal, currently only has effect on the BootstrapSample option. TRUE by default. |
| Category                     | major category to which the systematic belongs (instrumental, theory, ttbar, ...): used to split pulls plot for same category |
| SubCategory                  | minor category for the systematic, used to evaluate impact on POI per SubCategory in "i" step, defaults to Category setting if it is used, otherwise defaults to "Uncategorised", do not use "Gammas", "FullSyst", or "combine" as SubCategory names (reserved for special functionality) |
| CombineName                  | A unique string for each systematic that you want to combine into a single systematic (e.g.) envelope. This needs to be set for every systematic that needs to be combined. The code will then combine all systematic and modify the _first one_ and then set all other systematics manually to zero. This is executed during the `b` step. |
| CombineType                  | Can be `ENVELOPE` or `STANDARDDEVIATION`. Tells the code how to combine the systematics with the same `CombineName`. `STANDARDDEVIATION` does OneSided symmetrisation while for `ENVELOPE` it is recommended to use the symmetrisation for the individual components. |
| HistoPathUp                  | only for option HIST, for HISTO or SHAPE systematic: histogram file path for systematic up variation |
| HistoPathDown                | only for option HIST, for HISTO or SHAPE systematic: histogram file path for systematic down variation |
| HistoPathSufUp               | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic up variation |
| HistoPathSufDown             | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic down variation |
| HistoFile(s)Up               | only for option HIST, for HISTO or SHAPE systematic: histogram file name(s) for systematic up variation |
| HistoFile(s)Down             | only for option HIST, for HISTO or SHAPE systematic: histogram file name(s) for systematic down variation |
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
| Symmetrisation               | can be ONESIDED, TWOSIDED, ABSMEAN and MAXIMUM (...); for no symmetrisation, skip the line; ONESIDED = only one variation provided (e.g. up variation), down variation will be added as mirrored version, TWOSIDED = (up-down)/2 variation is calculated bin-by-bin, this is used as up variation and then mirrored to down variation, ABSMEAN = ((abs(up)+abs(down))/2) can be used when both variations have the same sign, USE WITH CAUTION, MAXIMUM = take variations with larger abs value w.r.t nominal and take this in each bin, can be used when both variations have the same sign, USE WITH CAUTION |
| Smoothing                    | smoothing code to apply; use 40 for default smoothing; for no smoothing, skip the line |
| OverallUp                    | for OVERALL systematic: the relative "up" shift (0.1 means +10%) |
| OverallDown                  | for OVERALL systematic: the relative "down" shift (-0.1 means -10%) |
| ScaleUp                      | for OVERALL, HISTO or SHAPE systematic: scale difference between "up" and nominal by a factor, or different factors for different regions (with the syntax `region1:1.2,region2:0.9`) |
| ScaleDown                    | for OVERALL, HISTO or SHAPE systematic: scale difference between "down" and nominal by a factor, or different factors for different regions (with the syntax `region1:1.2,region2:0.9`) |
| ReferenceSample              | if this is specified, the syst variation is evaluated w.r.t. this reference sample (often a GHOST sample) instead of the nominal, and then the relative difference is propagated to nominal; NOTE: also the overall relative difference is propagated |
| ReferenceSmoothing           | if this is specified, the syst variation is smoothed wrt a specified sample (Must appear in `Samples` for this syst) and then the smoothed variations are copied bin by bi to all other samples specified. Useful when `Morphing` is used |
| ReferencePruning             | if this is specified, the syst variation is pruned wrt a specified sample (Must appear in `Samples` for this syst) and then the same pruning is applied to all specified samples |
| DropShapeIn                  | specify regions or samples where you want the smoothing / pruning to be forced to drop the shape and keep only norm. When `all` is used, the shape is dropped for all regions and all samples |
| DropNorm                     | the same as the previous one, but to drop the norm and keep only the shape |
| DropNormSpecial              | Only effective when `CombineName` is used to combine special systematics. A list of region names can be specified to drop the normalisation for the combined systematic in those regions. This option has to be used for at least one of the systematics that are being combined |
| KeepNormForSamples           | list of samples (or sum of samples, in the form smp1+smp2), comma separated, for which the systematic gets shape only in each region |
| DummyForSamples              | list of samples, comma separated, for which the systematic gets created as a dummy one (syst variation = nominal); useful when used in combination with KeepNormForSamples |
| NoPruning                    | Force the code to skip pruning for the systematic |
| PreSmoothing                 | if set to TRUE, a TH1::Smooth-based smoothing is applied, prior to the usual smoothing (if set) |
| SmoothingOption              | if not set smoothing option from Job option is taken. Choose which smoothing option to use for this systematic (this will overwrite the general smoothing option for this), allowed parameters are: MAXVARIATION (default), TTBARRESONANCE (see also [FAQ section](#faq)), COMMONTOOLSMOOTHMONOTONIC, COMMONTOOLSMOOTHPARABOLIC, KERNELRATIOUNIFORM, KERNELDELTAGAUSS or KERNELRATIOGAUSS. |
| SubtractRefSampleVar         | if set to TRUE, the relative variation of the ReferenceSample will be linearly subtracted from the relative variation of each affected sample, for the same systematic - this is relevant e.g. for Full JER SmearingModel, where data would be the reference sample |
| HistoPathUpRefSample         | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file path for systematic up variation |
| HistoPathDownRefSample       | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file path for systematic down variation |
| HistoPathSufUpRefSample      | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram file names for systematic up variation |
| HistoPathSufDownRefSample    | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram file names for systematic down variation |
| HistoFile(s)UpRefSample      | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file name(s) for systematic up variation |
| HistoFile(s)DownRefSample    | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file name(s) for systematic down variation |
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
| Decorrelate                  | decorrelate systematic, can take values REGION (decorrelate across regions), SAMPLE (decorrelate across samples), SHAPEACC (decorrelate shape and acceptance effects); can be used to change behaviour of a systematic without having to re-run the n or b step |
| ForceShape                   | Can be: `NOSHAPE` (default) - no changes to the uncertainty, `LINEAR` will create a linear shape from left up variation to bottom down variation ased on the original up and down variations or `TRIANGULAR` that will create a triangular shape that is increasing until the middle of the distribution and then decreasing back to 0. The other variation is symmetrised. |



## Multi-fit
These options are for multi-fits, performed with action `m`.

### Multi-fit `Job` block settings
| **Option** | **Function** |
| ---------- | ------------ |
| Label            | the label which will be shown on plots |
| OutputDir        | the name of the output directory |
| LumiLabel        | the luminosity label that will be shown on the plots |
| CmeLabel         | the center of mass energy label that will be shown on the plots |
| CombiLabel         | the name of the combination that will be shown on the plots, default is "Combined" |
| SaveSuf          | added to file name of histograms, for usage with hupdate (equivalent to command line option) |
| ShowObserved     | can be TRUE or FALSE, flag to turn on/off the observed values on the plots |
| LimitTitle       | the title for limit that will be shown on the plots |
| POITitle         | the title of the POI that will be shown on X axis  |
| CompareLimits    | can be TRUE or FALSE, flag to compare to Limit values |
| ComparePOI       | can be TRUE or FALSE, flag to compare to POI values |
| ComparePulls     | can be TRUE or FALSE, flag to compare to pulls values |
| PlotCombCorrMatrix | can be set to TRUE or FALSE, flag to build correlation matrix from the combined systematics |
| CorrelationThreshold | Threshold used to draw the correlation matrix (only systematics with at least one correlation larger than than draw) (0.05:5%) |
| UseGammasForCorr | If set to `TRUE` will add gammas into correlation matrix plot. Default is `FALSE` |
| Combine          | can be TRUE or FALSE, set to TRUE if you want to perform actual combination (followed by `mwf`) |
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
| NuisParListFile  | Name of file containing list of nuisance parameters, with one parameter per line, and names just like in the `Fits/*txt` file. The order will be used for the plots created with `ComparePulls`. |
| PlotSoverB       | if set to TRUE will plot signal over background plots |
| SignalTitle      | a title of the signal for the plots |
| FitResultsFile   | a name of the file with fit results |
| LimitsFile       | a name of the file with limits results |
| BonlySuffix      | a suffix of the background only fits |
| ShowSystForPOI   | can be TRUE or FALSE, set to TRUE if you want to show systematics for POI |
| GetGoodnessOfFit | can be TRUE or FALSE, set to TRUE to get chi2/NDF for the fit. By default, the "SaturatedModel" is used (which is the only correct way how to estiamte the chi2/NDF). If you want to use the old test using Asimov, set GetGoodnessOfFit to TRUE and set SaturatedModel to FALSE |
| doLHscan         | comma separated list of NP(or POIs) to run LH scan, if first parameter is "all" it will be run for all NP |
| do2DLHscan       | produces 2D likelihood scan between the chosen parameters. Syntax: "paramX1,paramY1:param X2,paramY2". Warning takes long time. You can reduce the number of steps via `LHscanSteps`. Alternatively you can split up the 2D scan in slices with `Parallel2Dscan` |
| LHscanMin        | minimum value for the LH scan on x-axis (default is Norm min). This also effect the x-axis in a 2D scan |
| LHscanMax        | maximum value for the LH scan on x-axis (default is Norm max). This also effect the x-axis in a 2D scan |
| LHscanSteps      | number of steps on the LH scan (default is 30) . This also effect the x-axis in a 2D scan ||
| LHscanMinY       | minimum value for the 2D-LH scan on y-axis (default it Norm min) |
| LHscanMaxY       | maximum value for the 2D-LH scan on y-axis (default is Norm max) |
| LHscanStepsY     | number of steps on the LH scan in y-direction (default is 30, but if not specified it uses LHscanSteps) |
| PlotOptions      | same as for "standard" fits |
| Logo             | can be TRUE or FALSE, use TRUE to show TRExFitter logo |
| DebugLevel       | set level of debug output |
| RunROOTMacros    | can be TRUE or FALSE, set to TRUE to run the common scripts in root interpreter in stead of running the directly compiled version (FALSE, default) |
| POILabel         | name of the POI shown on plots, default is `#\mu` |
| POINominal       | value of the nominal (SM) prediction for POI, default is `1` |
| ShowTotalOnly    | If set to TRUE will show only total uncertainty on the POI plots. Default is FALSE |
| POIInitial       | Sets the initial value of the POI for the fit. Default is 1. |


### Multi-fit `Fit` block settings
| **Option** | **Function** |
| ---------- | ------------ |
| Options          | additional options, accepting only float as arguments - useful for adding your functionalities & flags in a quick way, since they need minimal changes in the code) ... |
| Label            | the label of the values from this config that will be shown on the plots |
| LoadSuf          |
| ConfigFile       | the path to the config file that you want to combine/compare |
| Workspace        | the path to the workspace that you want to combine/compare |
| ShowObserved     | can be TRUE or FALSE, set to TRUE to show the observed values of POI |
| FitResultsFile   | the path to the file with fit results |
| LimitsFile       | the path to the file with limits results |
| POIName          | the name of the POI |
| Directory        | the path to the directory |
| InputName        | the name of the input |
| UseInFit         | Set to `TRUE` by default, if set to `FALSE` will use the config in comparison plots, but will not be used to create the workspace and will not be used in the fit |
| UseInComparison         | Set to `TRUE` by default, if set to `FALSE` will not use the config in comparison plots. |



### Multi-fit `Limit` block settings
| **Option** | **Function** |
| ---------- | ------------ |
| LimitType                    | can be ASYMPTOTIC or TOYS (the latter is not yet supported) |
| LimitBlind                   | can be TRUE or FALSE (TRUE means that ALL regions are blinded) |
| POIAsimov                    | value of the POI to inject in the Asimov dataset in LimitBlind is set to TRUE |
| SignalInjection              | if set to TRUE, expected signal with signal injection is evaluated |
| SignalInjectionValue         | Value for the injected signal |
| ParamName                    | Name for the parameter in the output ROOT file |
| ParamValue                   | Value of the parameter in the output file (e.g. 172.5 for top mass) |
| OutputPrefixName             | Prefix for the output ROOT file |
| ConfidenceLevel              | Confidence level for the CLs. Default is 0.95 |


### Multi-fit `Significance` block settings
| **Option** | **Function** |
| ---------- | ------------ |
| SignificanceBlind            | can be TRUE or FALSE (TRUE means that ALL regions are blinded) |
| POIAsimov                    | value of the POI to inject in the Asimov dataset in LimitBlind is set to TRUE |
| ParamName                    | Name for the parameter in the output ROOT file |
| ParamValue                   | Value of the parameter in the output file (e.g. 172.5 for top mass) |
| OutputPrefixName             | Prefix for the output ROOT file |

### UnfoldingSample block
| **Option** | **Function** |
| ---------- | ------------ |
| Title                    | for plots |
| FillColor                | for plots |
| LineColor                | for plots |
| ResponseMatrixFile(s)        | same as for Job but per sample |
| ResponseMatrixName(s)        | same as for Job but per sample |
| ResponseMatrixPath(s)        | same as for Job but per sample |
| AcceptanceFile(s)            | same as for Job but per sample |
| AcceptanceName(s)            | same as for Job but per sample |
| AcceptancePath(s)            | same as for Job but per sample |
| SelectionEffFile(s)          | same as for Job but per sample |
| SelectionEffName(s)          | same as for Job but per sample |
| SelectionEffPath(s)          | same as for Job but per sample |
| MigrationEffFile(s)          | same as for Job but per sample |
| MigrationEffName(s)          | same as for Job but per sample |
| MigrationEffPath(s)          | same as for Job but per sample |
| Regions                  | Comma separated list of regions |

### UnfoldingSystematic block
| **Option** | **Function** |
| ---------- | ------------ |
| Samples                   | UnfoldingSamples affected by this uncertainty |
| Region                    | Regions affected by this uncertainty |
| NuisanceParameter         | Used to correlate systematics |
| Type                      | HISTO/OVERALL/SHAPE/STAT |
| Title                     | for plots |
| Category                  | same as for Systematic |
| Category                  | same as for Systematic |
| Symmetrisation            | same as for Systematic |
| SmoothingOption            | same as for Systematic |
| ResponseMatrixFile(s)Up        | self-explanatory |
| ResponseMatrixName(s)Up        | self-explanatory |
| ResponseMatrixPath(s)Up        | self-explanatory |
| ResponseMatrixFile(s)Down        | self-explanatory |
| ResponseMatrixName(s)Down        | self-explanatory |
| ResponseMatrixPath(s)Down        | self-explanatory |
| AcceptanceFile(s)Up        | self-explanatory |
| AcceptanceName(s)Up        | self-explanatory |
| AcceptancePath(s)Up        | self-explanatory |
| AcceptanceFile(s)Down        | self-explanatory |
| AcceptanceName(s)Down        | self-explanatory |
| AcceptancePath(s)Down        | self-explanatory |
| SelectionEffFile(s)Up        | self-explanatory |
| SelectionEffName(s)Up        | self-explanatory |
| SelectionEffPath(s)Up        | self-explanatory |
| SelectionEffFile(s)Down        | self-explanatory |
| SelectionEffName(s)Down        | self-explanatory |
| SelectionEffPath(s)Down        | self-explanatory |
| MigrationFile(s)Up        | self-explanatory |
| MigrationName(s)Up        | self-explanatory |
| MigrationPath(s)Up        | self-explanatory |
| MigrationFile(s)Down        | self-explanatory |
| MigrationName(s)Down        | self-explanatory |
| MigrationPath(s)Down        | self-explanatory |

### Unfolding block
| **Option** | **Function** |
| ---------- | ------------ |
| MatrixOrientation         | Can be TRUTHONHORIZONTAL/TRUTHONVERTICAL, self explanatory |
| TruthDistributionPath         | Globally set the path, will be overwritten by TruthSample |
| TruthDistributionFile         | Globally set the path, will be overwritten by TruthSample |
| TruthDistributionName         | Globally set the path, will be overwritten by TruthSample |
| NumberOfTruthBins         | Number of truth bins (can be only one truth distribution but can have more reco distributions - regions) |
| Tau         | Tikhonov regularization parameter, use as e.g. `1:2,3:1.7`, this will set the parameter to bins 1 and 3, with values 2 and 1.7 respectively; alternatively, a single value can be given, which will be used for all the unfolding norm-factors |
| TitleX         | Name of the title on unfolded plots |
| TitleY         | Name of the title on unfolded plots |
| RatioYmax         | Used for SetRangeUser |
| RatioYmin         | Used for SetRangeUser |
| LogX         | Set to TRUE to use logarithmic scale for unfolded plots |
| LogY         | Set to TRUE to use logarithmic scale for unfolded plots |
| TitleOffsetX         | scale of the original offset (default is 1) |
| TitleOffsetY         | scale of the original offset (default is 1) |
| UnfoldingResultMin   | Minimum for the normalisation factors for bins (Default is 0) |
| UnfoldingResultMax   | Maximum for the normalisation factors for bins (Default is 2) |
| NominalTruthSample   | Name of the TruthSample that should be used as nominal for folding. Has to be set! |
| MigrationTitleX   | name of the title on migration/response matrix plots |
| MigrationTitleY   | name of the title on migration/response matrix plots |
| MigrationLogX   | for migration/response plots |
| MigrationLogY   | for migration/response plots |
| MigrationTitleOffsetX   | for migration/response plots |
| MigrationTitleOffsetY   | for migration/response plots |
| MigrationZmin   | for migration plots |
| MigrationZmax   | for migration plots |
| ResponseZmax   | for response plots |
| ResponseZmin   | for response plots |
| PlotSystematicMigrations   | if set to TRUE will plot migration/response plots for all systematics |
| MigrationText   | if set to TRUE will show numbers in migration/response matrix plots |

### TruthSample block
| **Option** | **Function** |
| ---------- | ------------ |
| Title         | for plots |
| FillColor         | for plots |
| LineStyle         | for plots |
| TruthDistributionPath         | folder path for truth distributions |
| TruthDistributionFile         | file path for truth distributions |
| TruthDistributionName         | name of the histogram in the file for truth distributions |
| UseForPlotting                | Can be set to `TRUE` (default) of `FALSE`, tells the code if the given truth distribution will appear on the final plot with unfolded data |
