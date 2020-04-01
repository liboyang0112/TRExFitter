# Useful cross-checks of the results

The PL technique is a very powerful tool but it should be used with caution.
So how do you know that the result obtained can be trusted? If you see none of the problems described above, this is already a good sign that the fit is stable.
TRExFitter allows few more option for the validity of the fit:

- Since commit `5c38594cb6a9830fbb1fadb3d7eccb41eae66817` (13.4.2019) a [saturated model](https://its.cern.ch/jira/browse/TTHFITTER-220) can be used to cross-checks the goodness of the fit to data.
- Other very powerful tool is the usage of `Validation` regions - regions (and/or actually different variables in the same region) that are not used in the fit, but the fit results (correlations, pulls, constraints) are propagated to those regions/distributions.
  You can even calculate $\chi^2$ including bin-by-bin and NP correlations via `PlotOptions: Chi2` and `GetChi2: STAT+SYST`.
  The post-fit agreement in the validation regions demonstrates that the pulls/constraints/correlation obtained from the fit are not only specific for the phase space or distributions in the fit but represent rather general features.
