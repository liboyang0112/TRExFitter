# Goodness-of-fit calculation

## Saturated model

Whenever a fit procedure is used, it is important to check the Goodness-of-fit (GoF) status. GoF is a metric that gives probability how well the fit model can describe the observed data.
If the GoF gives very small probability, the model should be checked.

For a long time in TRExFitter, a very ad hoc GoF test was implemented that comapred the likelihood values for the fit to data and a fit to Asimov dataset. This is obviously not a proper test.
The proper test is provided by [saturated model](http://www.physics.ucla.edu/~cousins/stats/cousins_saturated.pdf). In this test, the likelihoods (likelihood ratios) are compared.
One is obtained by fitting the data with the nominal model and the other one by still fitting the real data, but with the _saturated model, a model that has enough freedom that it will fit the data perfectly, without nuisance parameter pulls.
In other words, the saturated model is a modified model that matches the data. The likelihood ratio follows the $\chi^2$ distribution asymptotically (Wilks theorem), and thus can be used as a standard GoF test.
To further compare it to the $\chi^2$ test, the saturated model represents $\chi^2 = 0$ (perfect agreement).
Or it can be viewed as the constant term that is removed from the full likelihood
$$
 -2\ln L(\mu) = \chi^2(\mu) + const
$$

### TRExFitter implemetation

The saturated model as GoF test has been implemented in TRExFitter since tag `TtHFitter-00-04-05`. Since tag `TRExFitter-00-04-08` it is the default option when `GetGoodnessOfFit` is set to `TRUE`.
Technically, it is implemented byt giving _shape factors_, "normalisatio nfactors per bin", to each bin to allow the model to fit the data perfectly without the need to pull any NP. The minimisation procedure is then run to with the value of the likelihood that is then used in the GoF calculation.

Let us try to use this option now.
We will use a configuration file used in our CI tests. First produce the inputs:

```bash
trex-fitter n test/configs/FitExampleNtuple.config
```

And now, run the fit (and also create the workspace first)

```bash
trex-fitter wf test/configs/FitExampleNtuple.config
```

We get the fitted values and everything looks fine. Now, modify the `Fit` block of the config file and add the following option: `GetGoodnessOfFit: TRUE` . SInce we are using a very recent verion of TRExFitter, this will by default use the saturated model as GoF test.
Now, run the workspace creation and the fit again

```bash
trex-fitter wf test/configs/FitExampleNtuple.config
```

You should not see that the fit is run _twice_. First time the standard fit is run, then the fit with the saturated model is run. You shoud see that in the second step no pulls are present.
Check the lines that print the likelihood values, e.g.:

```bash
=== INFO::FittingTool::FitPDF:    -> Reduced Final value of the NLL = -998511.89450029970612376928
=== INFO::FittingTool::FitPDF:    -> Final value of the NLL = 1488.105500
=== INFO::FittingTool::FitPDF:    -> Final value of offset = 1446.032751
=== INFO::FittingTool::FitPDF:    -> Final NLL - offset = 42.072749
```

```bash
=== INFO::FittingTool::FitPDF:    -> Reduced Final value of the NLL = -998516.03562775580212473869
=== INFO::FittingTool::FitPDF:    -> Final value of the NLL = 1483.964372
=== INFO::FittingTool::FitPDF:    -> Final value of offset = 1446.032751
=== INFO::FittingTool::FitPDF:    -> Final NLL - offset = 37.931622
```

And finally, the GoF value is printed

```bash
=== INFO::TRExFit::PerformFit: ----------------------- GOODNESS OF FIT EVALUATION -----------------------
=== INFO::TRExFit::PerformFit:   NLL0        = 1483.964372
=== INFO::TRExFit::PerformFit:   NLL         = 1488.105500
=== INFO::TRExFit::PerformFit:   ndof        = 10
=== INFO::TRExFit::PerformFit:   dNLL        = 4.141127
=== INFO::TRExFit::PerformFit:   2dNLL/nof   = 0.828225
=== INFO::TRExFit::PerformFit:   probability = 0.601288
=== INFO::TRExFit::PerformFit: ----------------------- -------------------------- -----------------------
```

!!! warning "Fit does not converge"
    Due to the implementation of the saturated model, it is possible that this fit fails even when the "normal" fit succeeds.

