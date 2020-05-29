# Mixed data and MC fits for "realistic" predictions

## Introduction

In analyses where a certain blinding procedure is applied,
the Asimov data-set is used to test the performance of the fit in terms of expected precision on the determination of the POI
and of the expected constraints on the nuisance parameters.
In the same way, expected exclusion limits and expected significance can be obtained.
However, often some of the analysis bins are not blinded, and the user can look at the data there and test the performance of the fit by including the data in these bins.
For simplicity, we refer to these unblinded bins as control regions (CRs) and to the blinded bins as signal regions (SRs).
With the purpose of testing the pulls and constraints on the nuisance parameters, it is trivial to perform a background-only fit to these CRs only, but this simple setup has a few drawbacks:
mainly no information on the expected sensitivity are obtained and only part of the data constraining power is inspected (i.e. nuisance parameter constraints could be stronger when fitting the SRs in addition to the CRs).
Therefore, one usually has to combine the information from the Asimov data-set fit (and/or limit and significance extraction based on Asimov data-set only)
and from such B-only fit to CRs.
However, there is a way to merge the two steps together in a consistent way, in order to obtain the most accurate predictions for the expected results by keeping the SRs completely blind to data and exploit all the available real-data information from the CRs.
Such an approach is referred to here as "mixed data and MC fit" setup, which is available when performing a fit, when extracting a limit and when getting a significance.

## Basic idea

The idea is to first perform a B-only fit to the CRs, using real data.
This preliminary fit provides a set of fitted values for all the nuisance parameters
(or at least for all the nuisance parameters that affects the background prediction in the CRs).
These fitted nuisance-parameter values are then used to construct a "modified Asimov data-set" in the SRs, where no real data can be used for the fit or limit/significance extraction.
Finally, the needed computation (fit, limit or significance) is performed considering all the CRs and the SRs usually included in the analysis,
but using a "mixed" data-set where real data is used in the CRs and this "modified Asimov data-set" is used in the SRs.
See this [reference](http://atlas-stats-doc-dev.web.cern.ch/atlas-stats-doc-dev/recommendations/rec_diagnostics_checks/#postfit-expected-results).

!!! tip "NB"
    The described procedure is NOT equivalent to first fitting the real data in the CRs and then using the fitted nuisance parameter values to construct a pseudo-Asimov data-set in ALL the regions: if this is done, the fitted values of the nuisance parameters in the final, combined fit will be in general different from those obtained in the preliminary fit.
    Typically looser constraints will be obtained, as the constrain terms in the likelihood will push them closer to their nominal values.

## The setup

In order to tell `TRExFitter` to work in a such a "mixed fit" setup,
the easiest way is to specify which data to user for each region using the `DataType` option.
For example:

```
Region: "mySR"
  DataType: DATA
  ...

Region: "myCR"
  DataType: ASIMOV
  ...
```

If a fit/limit/significance computation is made including regions with different `DataType`,
the fitter will automatically recognize that the used data is "mixed",
and the two-step procedure described above will be performed.

!!! tip "NB"
    In order for the procedure to work properly, no `NPValues` nor `NPValuesFromFitResults` should be set in the `Fit` block.
    If any of these is set, the Asimov data-set will be produced the specified nuisance parameters set to these values and no preliminary
    fit on the CRs will be performed.

### Alternative setups

Alternatively, to obtain the same result, other options could be used.

1. If the flag `BlindSRs` could be set to `TRUE` under `Job`, there is no need to specify `DataType` for each region, and instead all the regions flagged with `Type: SIGNAL` will be considered as if they had `DataType: ASIMOV`.
    Therefore an equivalent setup as the above will be:
    ```
    Job: "myJob"
      BlindSRs: TRUE
      ...

    Region: "mySR"
      Type: SIGNAL
      ...

    Region: "myCR"
      Type: CONTROL
      ...
    ```

2. The steps above could be run "manually", to have more control on the procedure.
    The user should first perform a fit to the CRs only (a background-ony fit typically) by hand, by either specifying the regions to fit from the     command line or by setting the fit type to `CRONLY` (if all and only the regions flagged as `CONTROL` need to be fitted).
    This fit will store the nuisance parameter values in the fit-result txt file, under the directory `Fit/`.
    This file can then be used to create the pseudo-Asimov data-set in the signal regions.
    To do it, the user should make sure the proper `DataType` is set for each of the CRs and SRs (or that `BlindSRs` is set)
    and then add the option `NPValuesFromFitResults` pointing to the fit-result file just produced.
    Finally, a normal fit can be run (and no preliminary fit will be performed, as `NPValuesFromFitResults` is set).
    Notice that this is equivalent to (but often more convenient than) specify all the needed nuisance-parameter values with the option `NPValues` under `Fit`.
