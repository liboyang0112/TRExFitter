# Limit calculation

In case of no excess in a search analysis, upper limits on the cross-section of certain signal processes are often set. 
In `TRExFitter`, there is the possibility to compute such exclusion limits in the case of a single POI controlling the normalization of a certain signal sample. 
According to ATLAS recommendations, this is done using the [CLs method](https://www.pp.rhul.ac.uk/~cowan/stat/cls/CLsInfo.pdf).
The exclusion limit calculation in `TRExFitter` is left to the `CommonStatTools` package, which provides classes and scripts to perform various operations on `RooFit` workspaces by means of profile-likelihood based computations.

According to [Statistics forum recommnedations](http://atlas-stats-doc-dev.web.cern.ch/atlas-stats-doc-dev/recommendations/rec_stat_computations/),
for the computation of the exclusion limit of a certain signal strength parameter $\mu$, 
a one-sided test-statistic is defined (it is basically just the negative-log-likelihood-ratio for values of $\mu$ larger than the measured one). 
Using this $\mu$-dependent test-statistic, the proper value of $\mu$ in order to set a 95% confidence-level one-sided interval is then extracted (using the CLs method). 

To obtain an exclusion limit within `TRExFitter`, one needs first to make sure the `Limit` block in the configuration file is filled.
Then the step `l` will compute the significance, provided that a workspace was already created (with the step `w`). 
The exclusion limit computation will consider as the parameter of interest the POI indicated in the `Job` block, and by default will take into account the real data (if present) both for the computation of the "observed" significance of the "expected" significance.
The output of the significance computation (performed internally by the `runAsymptoticsCLs` function in `CommonStatTools`) is stored in a `TTree` object saved in a ROOT file (`JOB-NAME/Limits/asymptotics/LIMIT-NAME_CL95.root`). 
Here one can see the values for the observed CLs limit (`obs_upperlimit`), the corresponding expected limit (`exp_upperlimit`),
plus a lot of more information.
Of particular interest are the values of the expected exclusion limit "$\pm$ 1 sigma" and "$\pm$ 1 sigma", which can be used to draw the green and yellow areas usually seen in the so-called "Brazilian plots". 
These are stored in the `exp_upperlimit_plus1`, `exp_upperlimit_plus2`, `exp_upperlimit_minus1` and `exp_upperlimit_minus2` branches. 

## Expected limits

Similarly to the case of the significance computation, the expected limits calculated by `TRExFitter` (or actually be `CommonStatTools`) are not "signal-blind". 
In order to build the Asimov data-set used to extract the expected exclusion limit, the values of the nuisance parameters fitted in the real data are used. 
This gives a more realistic value to compare with the observed one, especially in cases where large nuisance-parameter pulls change the background prediction significantly with respect to the pre-fit case. 
The user can force `TRExFitter` not to take into account the real data in this way and to stick to pre-fit predictions for the background, 
by simply setting `LimitBlind: TRUE`.

## Signal injection

By default, the POI is kept fixed to zero in the Asimov-data creation when computing the expected limit, 
so this expected limit is in fact "expected in the no-signal hypothesis". 
Optionally, another computation of the expected limit can be performed, by instead setting the value of the POI used in the Asimov-dataset creation to a non-zero value, by setting `SignalInjection: TRUE` and `SignalInjectionValue: SOME-NON-ZERO-NUMBER` in the `Limit` block. 
The resulting expected exclusion limit (only the central value, not the 1 and 2 sigma bands) is stored in the branch `inj_upperlimit` of the output `TTree`.
