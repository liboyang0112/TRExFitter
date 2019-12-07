# Job block settings

**Quick navigation**:
- Standard Fit
    - [Job block](Job.md)
    - [Fit block](Fit.md)
    - [Limit block](Limit.md)
    - [Significance block](Significance.md)
    - [Options block](Options.md)
    - [Region block](Region.md)
    - [Sample block](Sample.md)
    - [NormFactor block](NormFactor.md)
    - [ShapeFactor block](ShapeFactor.md)
    - [Systematic block](Systematic.md)
- Multi-Fit
    - [Job block](Job-multifit.md)
    - [Fit block](Fit-multifit.md)
    - [Limit block](Limit-multifit.md)
    - [Significance block](Significance-multifit.md)

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
| SystControlPlots             | if set to TRUE, plots showing the shape effect of a given systematic (before and after smoothing/symmetrisation) will be produced |
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
| SystErrorBars                | TRUE by default to add stat error bars to syst variations in syst plots, set to FALSE to disable |
| SummaryPlotRegions           | list of regions to be shown in summary plot (useful to specify a custom order) |
| FixNPforStatOnly             | if set to TRUE, when running stat-only (with either of the two options) also the norm factors other than the POI are kept fixed |
| InputFolder                  | specify it to read fit input histograms from a different directory than `<jobName>/Histograms/` |
| InputName                    | specify it to read fit input histograms from files with different name than `<jobName>_blabla.root` |
| OutputDir                    | specify it to write everything in a different directory than `<jobName>` |
| WorkspaceFileName            | if specified, an external ws can be used as input for fitting (not 100% supported) |
| KeepPruning                  | if set to TRUE, the first time the ws is created (option w) a Pruning.root file is created under `<jobName>/` and used for future operations to skip pruned systematics (makes operations much faster in case many syst are pruned) |
| AtlasLabel                   | to specify Internal, Preliminary, etc... If set to `none` the whole label will be removed |
| CleanTables                  | if set to TRUE, a cleaned version of the tex tables is created (basically removing the "#") - to be expanded |
| SystCategoryTables           | if set to TRUE, additional syst tables with systematics grouped by category are created |
| SummaryPlotYmax              | if set, it will force the summary plot to use this value as max y-maxis value |
| SummaryPlotYmin              | if set, it will force the summary plot to use this value as min y-maxis value |
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
| BootstrapSyst                | Specifies the systematics for which systematic should be the bootstrap done. |
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
