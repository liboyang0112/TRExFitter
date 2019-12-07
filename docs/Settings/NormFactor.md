# NormFactor block settings

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
| Tau                          | if set, it will add a constraint term for this norm factor (needed for unfolding regularization); the value is the value of the Tikhonov "tau" parameter: the largest the value, the stronger the regularization (smaller pre-fit uncertainty) |