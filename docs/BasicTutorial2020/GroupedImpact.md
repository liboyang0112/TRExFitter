# Grouping sources of uncertainty together

Besides the [ranking feature](Ranking.md), `TRExFitter` includes another way of calculating how much certain nuisance parameters "matter".
The feature discussed here is also called "grouped impact".
It is particularly useful to evaluate the uncertainty on a parameter of interest (POI) due to a group of nuisance parameters (NPs).

## Impact of (a group) of nuisance parameters

It works as follows.
The nominal `f` fit results in a POI uncertainty $\Delta \mu$.
When removing any nuisance parameters from the fit, the resulting new POI uncertainty $\Delta \mu^\prime$ should generally decrease: $\Delta \mu^\prime < \Delta \mu$.
If the uncertainty decreases a lot, then the nuisance parameter that was removed contributed a lot to the total uncertainty.
This procedure of removing parameters (or groups of parameters) can be used to define the impact of the parameter (or group) via quadrature subtraction of the POI uncertainties: $\sqrt{(\Delta \mu)^2 - (\Delta \mu^\prime)^2}$.
When splitting up all NPs into groups, the quadrature sum of the of the uncertainty contributions from all groups in general will be different from the total POI uncertainty.
This is due to correlations of NPs between the groups.

## How to perform a calculation

This method of evaluating the contribution to the total POI uncertainty from a group of nuisance parameters requires a few things.
First, the parameters need to be grouped together.
Two groups are defined by default: `Gammas` (MC statistical uncertainty) and `FullSyst` (full systematic impact, statistical component subtracted).
Additional groups can be defined in the configuration file, using the options `Category` and `SubCategory`.
If set, `SubCategory` overwrites the setting of `Category`.
A `NormFactor` defaults to `Category` called `NormFactors`.
If neither `Category` nor `SubCategory` are set, the group of a `Systematic` defaults to `Uncategorised`.
With all parameters grouped together as desired, the calculation of the impact of a group via quadrature subtraction can be performed.

The `TRExFitter` action `i` is used to perform the calculation.
It performs the nominal fit first, and then performs two more fit: a repeat of the nominal fit (see [Technical details](#technical-detail) below), and then a fit with the relevant group of nuisance parameters held constant (at their best-fit values).
Make sure to set `DebugLevel` to at least `1` to see relevant output.
Everything can be run in one step via

```bash
trex-fitter i <config>
```

If you have many groups, this might take a while and in that case it could be a good idea to parallelize:

```bash
trex-fitter i <config> GroupedImpact="Gammas"
```

If you parallelize, combine everything into a single file in the end via

```bash
trex-fitter i <config> GroupedImpact="combine"
```

The results are collected in `Fits/GroupedImpact.txt`.
Individual group results are calculated in a file with the name appended by the group name.
The file collecting all results contains one line per group, with lines looking like this:

```
FullSyst    9.20925  ( +8.75606, -10.4972 )
```

It lists the group name, followed by the result of the quadrature subtraction when using the symmetrized POI uncertainty, and in brackets are the results of the quadrature subtraction when using the uncertainty in up/down direction.
You might encounter `nan` results in the file.
Be on the lookout for errors related to this:

```
=== ERROR::FittingTool::FitExcludingGroup: uncertainty has increased for Gammas! please check the fit
```

See [Technical details](#technical-detail) below for more on this.

## Example

Let's try it out on the fit used in the [ranking example](Ranking.md).
Run

```bash
trex-fitter nwfr config/ttH2015_ljets.config
```

to set everything up, using the configuration file provided within the `TRExFitter` repository.
Next, we run the calculation all in one step:

```bash
trex-fitter i config/ttH2015_ljets.config
```

Set `DebugLevel: 1` to see the output, a part of it is replicated below:

```
=== INFO::FittingTool::GetGroupedImpact: -----------------------------------------------------
=== INFO::FittingTool::GetGroupedImpact: replicated nominal fit
=== INFO::FittingTool::GetGroupedImpact: POI is:   -4.073448 +/- 11.264900    ( +11.252870, -12.055149 )
=== INFO::FittingTool::GetGroupedImpact: -----------------------------------------------------
=== INFO::FittingTool::GetGroupedImpact: category: FullSyst (fixed to best-fit values for fit)
=== INFO::FittingTool::GetGroupedImpact: POI is:   -4.430454 +/- 6.487507    ( +7.068127, -5.927449 )
=== INFO::FittingTool::GetGroupedImpact:            --> impact: 9.209246    ( +8.756064, -10.497236 )
=== INFO::FittingTool::GetGroupedImpact: -----------------------------------------------------
```

The nominal POI uncertainty is $\pm 11.3$, and the uncertainty due to systematic uncertainties is $\pm 9.2$.
Systematic uncertainties are the dominant source of uncertainty in this example.
Look further through the results to see that the `ttbar-model` group has the largest impact.

## Technical details

The example above produces a few results that are `nan`.
This happens when the POI uncertainty in the fit with parameters removed is larger than the uncertainty with all parameters.
Unfortunately this can happen if the fits do not find the global likelihood maximum.
Most commonly this happens with groups that have a very small impact, and the two fits converge to slightly different local minima.
To avoid this issue, the `i` action first runs the nominal fit (to set up the environment), but then repeats it to allow the repeated fit and the fit with the group held fixed to be as similar as possible in their configuration.
Still, `nan` results can often not be avoided for groups with extremely small impact.
In such cases it might make sense to group such sources of uncertainty together with others, given that their impact is negligible anyway.

Since `i` performs the nominal fit twice, you might expect that the results are the exact same.
There are slight differences, with currently unknown origin (see also [the JIRA](https://its.cern.ch/jira/projects/TTHFITTER/issues/TTHFITTER-229)).
The differences should however not be larger than those from fluctuating initial settings (via `SetRandomInitialNPval`).
If you see large changes between the two nominal fits, please get in contact.
