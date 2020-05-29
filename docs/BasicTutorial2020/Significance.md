# Significance calculation

In case of a search, it is important to be able to assess the significance of a possible excess in data compatible with a certain signal hypothesis, 
as well as the expected signal significance. 
The significance calculataion in TRExFitter is left to the CommonStatTools package, which provides classes and scripts to perform various operations on RooFit workspaces by means of profile-likelihood based comoutations.

The significance of an excess is computed by defining the so-called p0-value. 
In general, the p0 is the probability of getting a worse compatibility between the data and the background-only model than the observed one. 
The compatibility between the data and a certain model (i.e. the prediction with a certain set of parameter values) is given by the test-statistics, 
which in the case of the profile-likelihood formalism, is the profile-likelihood ratio. 
The significance, often expressed in terms of number of "sigmas", is then just the translation of such p0 value into Gaussian standard deviations.

The exact procedure to get the significance in the profile-likelihood framework is explained here: http://atlas-stats-doc-dev.web.cern.ch/atlas-stats-doc-dev/recommendations/rec_stat_computations/ (see the "Discovery testing" paragraph).

To perform a significance computation within TRExFitter, one needs first to make sure the `Significance` block in the configuration file is filled.
Then the step `s` will compute the significance, provided that a workspace was already created (with the step `w`). 
The significance computation will consider as parameter of interest the POI indicated in the `Job` block, and by default will take into account the real data (if present) both for the computation of the "observed" significance of the "expected" significance. 
The output of the significance computation (performed internally by the `runSig` function in CommonStatTools) is stored in a TTree object saved in a root file (`JOB-NAME/Significance/asymptotics/SIG-NAME_p0.root`). 
Here one can see the values for the observed p0-value and significance (`obs_pval`, `obs_sig`) and the corresponding expected ones (`med_pval`, `med_sig`). 

## Expected significance

By default, the expected significance is computed by building an Asimov dataset with the values of the nuisance parameters obtained from a fit to the real data, 
and therefore is not "blind" to the data. 
In order to compute a "blind" expected significance computation, one needs to use the option `SignificanceBlind: TRUE`. 
Furthermore, the value of the POI used to build such an Asimov dataset for the expected significance computation can be set using the `POIAsimov` option in the `Significance` block.
