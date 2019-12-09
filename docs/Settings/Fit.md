# Fit block settings

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
| FitType                      | can be SPLUSB (default) or BONLY to fit under the s+b or the b-only hypothesis |
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
| TemplateInterpolationOption  | Option only for morphing, tells the code which interpolation between the templates is used. Three possible options are available: LINEAR(default)/SMOOTHLINEAR/SQUAREROOT. All of these options basically use linear interpolation but SMOOTHLINEAR approximates it by integral of hyperbolic tangent and SQUAREROOT approximates it by $`\sqrt{x^2+\epsilon}`$ to achieve smooth transitions (first derivative) between the templates |
| BlindedParameters            | A comma separated list of POI/NPs that will be written as a hexadecimal number so it is not easy to read to not accidentally unblind. When at least one parameter is set the console output of the minimization is removed.
| DoNonProfileFitSystThreshold | When performing a NonProfileFit, systematics are not added to total if smaller than this threshold |