# Treatment of interference effects in TRExFitter

Interference effects between different processes cannot be neglected in some analyses. 
If the signal process interferes strongly with the main background, they can significantly alter the signal shape compared to a simple Breit-Wigner peak.
This is illustrated in the following figure for the interference between the loop induced production of heavy neutral Higgs bosons $`A/H`$ decaying to a top-antitop-quark pair
($gg \rightarrow A/H \rightarrow t\bar{t}$) and the background from SM $t\bar{t}$ production.

<img src="https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/EXOT-2016-04/fig_01.png"
     alt="Interference"
     width="400"
     />

## Extension of the signal model

For a signal hypothesis test with interference effects taken into account, the signal model in the likelihood ($\mu \cdot S + B$ without interference) can be extended as follows:

```math
\mu \cdot S + \sqrt{\mu} \cdot I + B = (\mu - \sqrt{\mu}) \cdot S + \sqrt{\mu} \cdot (S+I) + B
```
The left-hand (right-hand) side of the equation is used if pure interference $I$ (signal-plus-interference $S+I$) histogram are available in the template fit.

The $\mu$ dependence of the extended signal model can be easily implemented in the `TRExFitter` config file 
using the [expression syntax](../AdvancedTutorial2020/Expression.md). Thus, the right-hand side of the above equation could be expressed as follows in the config file,
in which sqrt_mu is taken to be the parameter of interest (POI):

```yaml
Job: "OffsetMethodForInterference"
  ...
  POI: "sqrt_mu"
  ...
 
NormFactor: "sqrt_mu"
  Samples: "SplusI"
  Nominal: 1
  Min: 0
  Max: 5

NormFactor: "norm_fact_pure_signal"
  Samples: "PureS"
  Expression: (sqrt_mu^2 - sqrt_mu):sqrt_mu[1,0,5]
```



## Treatment of histograms with negative bin content

Bins with negative content, which appear often in histograms with interference (both $S+I$ or $I$) cannot easily be handled by the underlying statistical tools (in particular `HistFactory`).
`TRExFitter` contains various layers of protection against negative bin content, which it is usually considered unphysical.

Note that, in principle, the above signal model could be rewritten in terms of the inclusive process histogram ($S+I+B$), which has only positive-entry bins. 
However, the generation of the inclusive $S+I+B$ process for each signal parameter point is usually prohibitively CPU-expensive (due to the background cross-section being much larger
than the signal cross-section), hence it is usually preferable, from the MC-generation point of view, to work with either $S+I$ or $I$ histograms.

### Offset method

A simple way to process histograms with negative bin content in the context of the extended signal model described above is to simply add an offset histogram 
to the respective histograms ($I$ or $S+I$) such that the resulting histogram has only positve entries.

To compensate for the offset, an equivalent "counterterm" histogram
has to be added to the template fit. This template histogram is identical to the offset histogram but is scaled with a negative normalisation factor.

In short, we replace the term $\sqrt{\mu} \cdot (S+I)$ in the signal model by

```math
\sqrt{\mu} \cdot [(S+I)+\texttt{offset}] - \sqrt{\mu} \cdot \texttt{counterterm}
```

where $(S+I)$+\texttt{offset} now has only positive entries and \texttt{counterterm} is identical to the offset histogram.
The same procedure can be applied to the pure interference histogram $I$ (left-hand side of the extended signal model equation).

The offset histogram can either be an artifical flat offset histogram (with bin contents chosen appropriately to compensate the negative content in the $S+I$ or $I$ histograms).
Alternatively, one of the background histograms, for example the one for the interfering background, can be chosen as an offset (and counterterm) histogram.

### Implementation of offset method in the config file

In the `TRExFitter` config file, the offset histogram is added as a `GHOST` sample, i.e. it has no associated statistical or systematic uncertainties.

```yaml
Sample: "SplusI_Offset"
  Type: GHOST
  HistoFile: "bkg1"
...

Sample: "SplusI"
  Type: SIGNAL
  Title: "Signal+Interference"
  FillColor: 632
  LineColor: 632
  HistoFile: "siginter"  
  SeparateGammas: True
  AddSample: "SplusI_Offset"
```

The counterterm is added as an extra sample with MC statistical errors not considered (see below).

```yaml
Sample: "SplusI_CounterTerm"
  Type: SIGNAL
  Title: "Counterterm"
  HistoFile: "bkg1"
  UseMCstat: False 
```

Its normalisation factor has the same magnitude as that for the respective $S+I$ or $I$ sample but with a negative sign:

```yaml
NormFactor: "minus_sqrt_mu"
  Samples: "SplusI_CounterTerm"
  Expression: -sqrt_mu:sqrt_mu[1,0,5]
```

An example config file with a minimal working example implementing the offset method can be found here:
[test/configs/FitExampleInterferenceOffset.config](https://gitlab.cern.ch/TRExStats/TRExFitter/-/tree/master/test/configs)

#### MC statistical errors with the offset method

The offset approach involves different MC statistical errors compared to the equivalent implementation without the offset histograms and hence different pulls and constraints
for the `gamma` nuisance parameters. For analyses that are completely systematics dominated, the resulting changes may be negligible.
However, it is recommended to correct the MC statistical uncertainty for the $S+I$ or $I$ histograms by hand. This is achieved as follows:

* Disable the MC statistical errors for the offset and counterterm histograms (`UseMCStats: False`). 
* Additionally, add separate `gamma` factors for the $S+I$ or $I$ histograms (`SeparateGammas: True`).

While this approach doubles the number of `gamma` nuisance parameters for an analysis, this approach is actually the more correct approach compared to the standard approach of having
only one `gamma` factor per bin (combining all samples), see the original [Barlow-Beeston method](https://www.sciencedirect.com/science/article/pii/001046559390005W).

#### Smoothing and pruning with the offset method

The addition of an offset histogram can affect the results of the smoothing and pruning procedures applied to the systematics histograms before the template fit.
The effect is expected to be small for most applications, however, analysers should be careful to make sure that this is indeed the case.
The cross-checks in the following sections are therefore recommended for interference analyses using the offset method.

## Recommended cross-checks for analysers

To ensure that the offset method is implemented correctly and smoothing/pruning results are unchanged after the addition of the offset histograms, the following cross-checks 
should be performed by every interference analysis using the offset method.

* Perform a simple `S+B` fit (dropping the `S+I` or `I` component from the signal model).
* Compare the results of this fit (best-fit values, ranking, impact) with those from an equivalent fit in which the offset is added to the signal (`S`) component.
* Significant differences could point to differences in smoothing and/or pruning.
* This test should be performed for (at least) 1-2 representative signal points (e.g. one high-mass and one low-mass point).
