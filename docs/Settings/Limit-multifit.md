# Multi-Fit Limit block settings

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
| LimitType                    | can be ASYMPTOTIC or TOYS (the latter is not yet supported) |
| LimitBlind                   | can be TRUE or FALSE (TRUE means that ALL regions are blinded) |
| POIAsimov                    | value of the POI to inject in the Asimov dataset in LimitBlind is set to TRUE |
| SignalInjection              | if set to TRUE, expected signal with signal injection is evaluated |
| SignalInjectionValue         | Value for the injected signal |
| ParamName                    | Name for the parameter in the output ROOT file |
| ParamValue                   | Value of the parameter in the output file (e.g. 172.5 for top mass) |
| OutputPrefixName             | Prefix for the output ROOT file |
| ConfidenceLevel              | Confidence level for the CLs. Default is 0.95 |