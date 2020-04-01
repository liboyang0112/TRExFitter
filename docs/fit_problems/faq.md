# FAQ

???+ question "Can we insert correlations for NPs?"
    No.
    this is not possible as the fit finds the correlations as described in [Minimisation procedure](technical_aspects.md#minimisation-procedure).
    Also this is problematic due to the split of the effects.
    All we can do it apply 100% correlation (fully correlating two NPs) - this can be done by providing same `NuisanceParameter` string for NPs that you want to fully correlate

???+ question "Can we insert correlations for NFs?"
    This is possible!
    There is a difference between correlating NPs and NFs as NFs do not have a prior term.
    Have a look at `Expression` option in the [TRExFitter readme](https://gitlab.cern.ch/TRExStats/TRExFitter).
    This can also be used to calculate some other properties and not only signal strengths (like ratio of signal strengths, fitting fractions directly, etc.)

???+ question "Why do I get different results every time I run identical setup?"
    This should not happen.
    It can happen only when you use `SetRandomInitialNPval` as the fit will start from a different point each time and thus the result can vary slightly (rule of thumb - the difference should not be larger than 0.01 on the POI).

???+ question "How should I fix empty bins?"
    TRExFitter fixes empty bins for you by adding some arbitrary (small) values, otherwise the fit would not work at all.
    *But this is not a good solution!* You should find a proper solution - merging small backgrounds, rebinning to get rid of zero/negative bin contents.

???+ question "What should I do with systematic uncertainties that have both up and down shift in the same direction and I want to symmetrise them?"
    TRExFitter will print warnings when this happens.
    Problem is that the symmetrisation procedure does not make sense when this happen.
    We provide two alternative symmetrisation options: `Symmetrisation: ABSMEAN` or `Symmetrisation: MAXIMUM`.
    However, there is no technical way how to fix this problem and these options are not enough.
    So you need to adjust the inputs, unfortunately.