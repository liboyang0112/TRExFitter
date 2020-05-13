# Using the `Expression` config option

## Introduction

It is technically not possible to add correlations of the NPs to the likelihood by hand. Correlations of the NPs are estimated from the fit during the likelihood minimisation step.
However, it is possible to correlate the normalisation factors with other normalisation factors.
This allows to modify the standard fits to non-standard ones, e.g. fraction fitting where the sum of some normalisation has to add up to one.
These correlations can be set via the `Expression` option.

!!! hint "Expression and morphing"
    The underlying functionality, the `AddPreprocess()` RooFit function, is used for the `Expression` part as well as for the fits with multiple templates.

## Example
