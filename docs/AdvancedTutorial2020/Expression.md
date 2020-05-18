# Using the `Expression` config option

## Introduction

It is technically not possible to add correlations of the NPs to the likelihood by hand.
Correlations of the NPs are estimated from the fit during the likelihood minimisation step.
However, it is possible to correlate the normalisation factors with other normalisation factors.
This allows to modify the standard fits to non-standard ones, e.g. fraction fitting where the sum of some normalisation has to add up to one.
These correlations can be set via the `Expression` option.

!!! hint "Expression and morphing"
    The underlying functionality, the `AddPreprocess()` RooFit function, is used for the `Expression` part as well as for the fits with multiple templates.

## Example

As an example use case for the `Expression` option is the W helicity measurement.
This measurement fits the fractions of the individual helicity templates (pure left-handed, pure longitudinal and pure right-handed), however, the fractions are not independent.
The fractions have to satisfy 

$$
F_{L} + F_{0} + F_{R} = 1
$$

for all allowed values of the fractions.
This leaves only two independent fractions.
The fit setup then needs three NormFactors, one for each fraction, but one of the NormFactors needs to be fixed from the other ones.

Have a look at the config file in `test/configs/FitExampleExpression.config` and focus on the line withe the Expression in the `NormFactor: "norm_left"` block

```bash
 Expression: (1.-norm_long-norm_right):norm_long[0.687,0,1],norm_right[0.002,0,1]
```

In this line, we tell the code to replace the `norm_left` normalisation parameter with `1 - norm_long - norm_right` and then set the intial and minimum/maximum values for the normalisation factors.

!!! hint "Note"
    The formula part of the expression can be anything that can be read by `TFormula`.

Now produce the histograms and run the fit

```bash
trex-fitter hwf test/configs/FitExampleExpression.config
```

You will see that the `norm_left` dissapeared from the results completely, but that is expected since this parameter no longer exists in the likelihood.

!!! hint "Flexibility"
    The `Expression` functionality provides high level of flexibility for different kind of measurements. Do not be afraid to experiment with it.
