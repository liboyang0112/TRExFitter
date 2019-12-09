# Region block settings

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
| NtuplePathSuff(s)            | only for option NTUP, the path suffix (or suffixes, comma-separated) where to find the ntuple files for this region |
| MCweight                     | only for option NTUP, the additional weight used in this region (for MC samples only) |
| Rebin                        | if specified, the histograms will be rebinned merging N bins together, where N is the argument (int) |
| Binning                      | if specified, the histograms will be binned according to the new binning specified, in the form like (0,10,20,50,100). If option AutoBin is set, use algorithms/functions or define the binning. Example - Binning: "AutoBin","TransfoD",5.,6. (TransfoF also available, 5. and 6. are parameters of the transformation). If used in background region and zSig!=0 (first parameter, =0 gives flat background) then need a comma separated list of backgrounds to use instead of signal to compute the binning. |
| Rebinning                    | if specified, the histograms will be rebinned according to the new binning specified, in the form like (0,10,20,50,100). Differently from the BInning option, this one performs the rebinning aftre the original histograms are created. This means that this option can changed (or removed) before running the b step. |
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
| XaxisRange                   | Manually call 'SetRangeUser()' on X axis. Needs two parameters(floats): min,max |
