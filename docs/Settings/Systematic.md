# Systematic block settings

**Quick navigation**:
- Standard Fit
    - [Job block](Job.md)
    - [Fit block](Fit.md)
    - [Limit block](Limit.md)
    - [Significance block](Significance.md)
    - [Options block](Options.md)
    - [Region block](Region.md)
    - [Sample block](Sample.md)
    - [NormFactor block](NormFactor.md)
    - [ShapeFactor block](ShapeFactor.md)
    - [Systematic block](Systematic.md)
- Multi-Fit
    - [Job block](Job-multifit.md)
    - [Fit block](Fit-multifit.md)
    - [Limit block](Limit-multifit.md)
    - [Significance block](Significance-multifit.md)

| **Option** | **Function** |
| ---------- | ------------ |
| Samples                      | comma-separated list of samples on which to apply the systematic |
| Regions                      | comma-separated list of regions where to apply the systematic |
| Exclude                      | comma-separated list of samples/regions to exclude |
| ExcludeRegionSample          | comma-separated list of region:sample to exclude |
| Type                         | can be HISTO, OVERALL, SHAPE (this refers to the HistFactory Shape Systematic, i.e. uncorrelated bin-by-bin) or STAT (this refers to auto-creation of one systematic from stat uncertainty for each bin of corresponding region - DEPRECATED) |
| Title                        | title of the systematic (will be shown in plots) |
| StoredName                   | if specified, will be used to read and write histograms in the root files under Histograms/ instead of the syst name; useful to decorrelate without re-creating histograms |
| NuisanceParameter            | if specified, this will be given to RooStats instead of the syst name; useful (and recommended) way to correlate systematics |
| IsFreeParameter              | if set to TRUE, the constraint will be a flat one instead of Gaussian (use with caution) |
| Category                     | major category to which the systematic belongs (instrumental, theory, ttbar, ...): used to split pulls plot for same category |
| SubCategory                  | minor category for the systematic, used to evaluate impact on POI per SubCategory in "i" step, defaults to Category setting if it is used, otherwise defaults to "Uncategorised", do not use "Gammas", "FullSyst", or "combine" as SubCategory names (reserved for special functionality) |
| CombineName                  | A unique string for each systematic that you want to combine into a single systematic (e.g.) envelope. This needs to be set for every systematic that needs to be combined. The code will then combine all systematic and modify the _first one_ and then set all other systematics manually to zero. This is executed during the `b` step. |
| CombineType                  | Can be `ENVELOPE` or `STANDARDDEVIATION`. Tells the code how to combine the systematics with the same `CombineName`. `STANDARDDEVIATION` does OneSided symmetrisation while for `ENVELOPE` it is recommended to use the symmetrisation for the individual components. |
| HistoPathUp                  | only for option HIST, for HISTO or SHAPE systematic: histogram file path for systematic up variation |
| HistoPathDown                | only for option HIST, for HISTO or SHAPE systematic: histogram file path for systematic down variation |
| HistoPathSufUp               | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic up variation |
| HistoPathSufDown             | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic down variation |
| HistoFile(s)Up               | only for option HIST, for HISTO or SHAPE systematic: histogram file name(s) for systematic up variation |
| HistoFile(s)Down             | only for option HIST, for HISTO or SHAPE systematic: histogram file name(s) for systematic down variation |
| HistoFileSufUp               | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic up variation |
| HistoFileSufDown             | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram file names for systematic down variation |
| HistoNameUp                  | only for option HIST, for HISTO or SHAPE systematic: histogram name for systematic up variation |
| HistoNameDown                | only for option HIST, for HISTO or SHAPE systematic: histogram name for systematic down variation |
| HistoNameSufUp               | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram names for systematic up variation |
| HistoNameSufDown             | only for option HIST, for HISTO or SHAPE systematic: suffix of the histogram names for systematic down variation |
| NtuplePath(s)Up              | only for option NTUP, for HISTO or SHAPE systematic: ntuple file path(s) for systematic up variation |
| NtuplePath(s)Down            | only for option NTUP, for HISTO or SHAPE systematic: ntuple file path(s) for systematic down variation |
| NtuplePathSufUp              | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file paths for systematic up variation |
| NtuplePathSufDown            | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file paths for systematic down variation |
| NtupleFile(s)Up              | only for option NTUP, for HISTO or SHAPE systematic: ntuple file name(s) for systematic up variation |
| NtupleFile(s)Down            | only for option NTUP, for HISTO or SHAPE systematic: ntuple file name(s) for systematic down variation |
| NtupleFileSufUp              | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file names for systematic up variation |
| NtupleFileSufDown            | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple file names for systematic down variation |
| NtupleName(s)Up              | only for option NTUP, for HISTO or SHAPE systematic: ntuple name(s) for systematic up variation |
| NtupleName(s)Down            | only for option NTUP, for HISTO or SHAPE systematic: ntuple name(s) for systematic down variation |
| NtupleNameSufUp              | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple names for systematic up variation |
| NtupleNameSufDown            | only for option NTUP, for HISTO or SHAPE systematic: suffix of the ntuple names for systematic down variation |
| SampleUp                     | if set, the syst variation will be built comparing the sample with another sample after all corrections are done; NB: can be used only if the syst affects one sample only |
| SampleDown                   | if set, the syst variation will be built comparing the sample with another sample after all corrections are done; NB: can be used only if the syst affects one sample only |
| WeightUp                     | only for option NTUP, for HISTO or SHAPE systematic: weight for systematic up variation; it will be multiplied with the `MCweight` defined in the Job and Region blocks, however it will not be multiplied with the `MCweight` defined in the sample block of the reference sample (for this, use `WeightSufUp` instead) |
| WeightDown                   | only for option NTUP, for HISTO or SHAPE systematic: weight for systematic down variation; it will be multiplied with the `MCweight` defined in the Job and Region blocks, however it will not be multiplied with the `MCweight` defined in the sample block of the reference sample (for this, use `WeightSufDown` instead) |
| WeightSufUp                  | only for option NTUP, for HISTO or SHAPE systematic: additional weight for systematic up variation (also multiplied with the `MCweight` acting on the nominal sample) |
| WeightSufDown                | only for option NTUP, for HISTO or SHAPE systematic: additional weight for systematic down variation (also multiplied with the `MCweight` acting on the nominal sample) |
| IgnoreWeight                 | only for option NTUP: if set, the corresponding weight (present in Job, Sample or Region) will be ignored for this systematic |
| Symmetrisation               | can be ONESIDED, TWOSIDED, ABSMEAN and MAXIMUM (...); for no symmetrisation, skip the line; ONESIDED = only one variation provided (e.g. up variation), down variation will be added as mirrored version, TWOSIDED = (up-down)/2 variation is calculated bin-by-bin, this is used as up variation and then mirrored to down variation, ABSMEAN = ((abs(up)+abs(down))/2) can be used when both variations have the same sign, USE WITH CAUTION, MAXIMUM = take variations with larger abs value w.r.t nominal and take this in each bin, can be used when both variations have the same sign, USE WITH CAUTION |
| Smoothing                    | smoothing code to apply; use 40 for default smoothing; for no smoothing, skip the line |
| OverallUp                    | for OVERALL systematic: the relative "up" shift (0.1 means +10%) |
| OverallDown                  | for OVERALL systematic: the relative "down" shift (-0.1 means -10%) |
| ScaleUp                      | for OVERALL, HISTO or SHAPE systematic: scale difference between "up" and nominal by a factor, or different factors for different regions (with the syntax `region1:1.2,region2:0.9`) |
| ScaleDown                    | for OVERALL, HISTO or SHAPE systematic: scale difference between "down" and nominal by a factor, or different factors for different regions (with the syntax `region1:1.2,region2:0.9`) |
| ReferenceSample              | if this is specified, the syst variation is evaluated w.r.t. this reference sample (often a GHOST sample) instead of the nominal, and then the relative difference is propagated to nominal; NOTE: also the overall relative difference is propagated |
| ReferenceSmoothing           | if this is specified, the syst variation is smoothed wrt a specified sample (Must appear in `Samples` for this syst) and then the smoothed variations are copied bin by bi to all other samples specified. Useful when `Morphing` is used |
| ReferencePruning             | if this is specified, the syst variation is pruned wrt a specified sample (Must appear in `Samples` for this syst) and then the same pruning is applied to all specified samples |
| DropShapeIn                  | specify regions or samples where you want the smoothing / pruning to be forced to drop the shape and keep only norm. When `all` is used, the shape is dropped for all regions and all samples |
| DropNorm                     | the same as the previous one, but to drop the norm and keep only the shape |
| KeepNormForSamples           | list of samples (or sum of samples, in the form smp1+smp2), comma separated, for which the systematic gets shape only in each region |
| DummyForSamples              | list of samples, comma separated, for which the systematic gets created as a dummy one (syst variation = nominal); useful when used in combination with KeepNormForSamples |
| PreSmoothing                 | if set to TRUE, a TH1::Smooth-based smoothing is applied, prior to the usual smoothing (if set) |
| SmoothingOption              | if not set smoothing option from Job option is taken. Choose which smoothing option to use for this systematic (this will overwrite the general smoothing option for this), allowed parameters are: MAXVARIATION (default), TTBARRESONANCE (see also [FAQ section](#faq)), COMMONTOOLSMOOTHMONOTONIC, COMMONTOOLSMOOTHPARABOLIC, KERNELRATIOUNIFORM, KERNELDELTAGAUSS or KERNELRATIOGAUSS. |
| SubtractRefSampleVar         | if set to TRUE, the relative variation of the ReferenceSample will be linearly subtracted from the relative variation of each affected sample, for the same systematic - this is relevant e.g. for Full JER SmearingModel, where data would be the reference sample |
| HistoPathUpRefSample         | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file path for systematic up variation |
| HistoPathDownRefSample       | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file path for systematic down variation |
| HistoPathSufUpRefSample      | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram file names for systematic up variation |
| HistoPathSufDownRefSample    | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram file names for systematic down variation |
| HistoFile(s)UpRefSample      | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file name(s) for systematic up variation |
| HistoFile(s)DownRefSample    | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram file name(s) for systematic down variation |
| HistoFileSufUpRefSample      | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram file names for systematic up variation |
| HistoFileSufDownRefSample    | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram file names for systematic down variation |
| HistoNameUpRefSample         | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram name for systematic up variation |
| HistoNameDownRefSample       | only for option HIST, for HISTO or SHAPE systematic: reference sample histogram name for systematic down variation |
| HistoNameSufUpRefSample      | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram names for systematic up variation |
| HistoNameSufDownRefSample    | only for option HIST, for HISTO or SHAPE systematic: reference sample suffix of the histogram names for systematic down variation |
| NtuplePath(s)UpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple file path(s) for systematic up variation |
| NtuplePath(s)DownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple file path(s) for systematic down variation |
| NtuplePathSufUpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple file paths for systematic up variation |
| NtuplePathSufDownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple file paths for systematic down variation |
| NtupleFile(s)UpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple file name(s) for systematic up variation |
| NtupleFile(s)DownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple file name(s) for systematic down variation |
| NtupleFileSufUpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple file names for systematic up variation |
| NtupleFileSufDownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple file names for systematic down variation |
| NtupleName(s)UpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple name(s) for systematic up variation |
| NtupleName(s)DownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample ntuple name(s) for systematic down variation |
| NtupleNameSufUpRefSample     | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple names for systematic up variation |
| NtupleNameSufDownRefSample   | only for option NTUP, for HISTO or SHAPE systematic: reference sample suffix of the ntuple names for systematic down variation |
| Decorrelate                  | decorrelate systematic, can take values REGION (decorrelate across regions), SAMPLE (decorrelate across samples), SHAPEACC (decorrelate shape and acceptance effects); can be used to change behaviour of a systematic without having to re-run the n or b step |