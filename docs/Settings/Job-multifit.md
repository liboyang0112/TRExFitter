# Multi-Fit Job block settings

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
| Label            | the label which will be shown on plots |
| OutputDir        | the name of the output directory |
| LumiLabel        | the luminosity label that will be shown on the plots |
| CmeLabel         | the center of mass energy label that will be shown on the plots |
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
| GetGoodnessOfFit | can be TRUE or FALSE, set to TRUE to get chi2/NDF for the fit |
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