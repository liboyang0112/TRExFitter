# Sample block settings

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
| Type                         | can be SIGNAL, BACKGROUND, DATA or GHOST; default is BACKGROUND; GHOST means: no syst, not drawn, not propagated to workspace |
| Title                        | title shown on the legends |
| TexTitle                     | title shown on tex tables |
| Group                        | if specified, sample will be grouped with other samples with same group and this label will be used in plots |
| HistoName(s)                 | valid only for option HIST; name(s) of histogram to read |
| HistoFile(s)                 | valid only for option HIST; which root file(s) to read (excluding the suffix ".root"); this will be combined with Fit->HistoPath to build the full path |
| HistoPath(s)                 | valid only for option HIST; it's the path(s) where the input root files containing the histograms are stored |
| HistoNameSuff(s)             | valid only for option HIST; suffix(es) for the name of histogram(s) to read |
| HistoFileSuff(s)             | valid only for option HIST; suffix(es) for the file name of histogram(s) to read |
| HistoPathSuff(s)             | valid only for option HIST; suffix(es) for the path of histogram(s) to read |
| NtupleName(s)                | valid only for option NTUP; name(s) of tree to read |
| NtupleFile(s)                | valid only for option NTUP; it's the file name(s) where the input ntuples are stored |
| NtuplePath(s)                | valid only for option NTUP; it's the path(s) where the input root files containing the ntuples are stored |
| NtupleNameSuff(s)            | valid only for option NTUP; suffix(es) for the name of tree to read |
| NtupleFileSuff(s)            | valid only for option NTUP; suffix(es) for the file name(s) of tree to read |
| NtuplePathSuff(s)            | valid only for option NTUP; suffix(es) for the path(s) of tree to read |
| FillColor                    | histogram fill color (not valid for data) |
| FillColorRGB                 | histogram fill color in RGB (not valid for data). This expects a triplet of RGB values between 0 and 255, e.g. `255,0,0`. If set, the FillColor option is ignored. |
| LineColor                    | histogram line color |
| LineColorRGB                 | histogram line color in RGB. This expects a triplet of RGB values between 0 and 255, e.g. `255,0,0`. If set, the LineColor option is ignored. |
| NormFactor                   | NormalisationFactor (free parameter in the fit); in the format \<name\>,nominal,min,max |
| ShapeFactor                  | ShapeFactor added |
| NormalizedByTheory           | set it to FALSE for data-driven backgrounds (MCweight, Lumi and LumiScale from Job and Region will be ignored) |
| MCweight                     | only for option NTUP, the additional weight used in this sample |
| Selection                    | valid only for option NTUP; additional selection for this region |
| Regions                      | set this to have the sample only in some regions |
| Exclude                      | set this to exclude the sample in some regions |
| LumiScale(s)                 | set this to scale the sample by a number; if more numbers are set, use a different one for each file / name / path... |
| IgnoreSelection              | if set, selection from Job and Region will be ignored |
| UseMCstat                    | if set to FALSE, makes the fitter ignore the stat uncertainty for this sample |
| UseSystematics               | has to be set to TRUE to allow systematics on the GHOST samples |
| MultiplyBy                   | if specified, each sample hist is multiplied bin-by-bin by another sample hist, in each of the regions |
| DivideBy                     | if specified, each sample hist is divided bin-by-bin by another sample hist, in each of the regions |
| AddSample(s)                 | if specified, each sample hist gets added bin-by-bin another sample hist, in each of the regions |
| SubtractSample(s)            | if specified, each sample hist gets subtracted bin-by-bin another sample hist, in each of the regions |
| Smooth                       | if set to TRUE, the nominal histograms are smoothed (based on TH1::Smooth but taking into account the original stat uncertainty) |
| AsimovReplacementFor         | only for GHOST samples; if set, the creation of custom Asimov data-set(s) is triggered; use as 'AsimovReplacementFor: "dataset","sample"', where "dataset" is the name of a custom Asimov dataset one wants to create (the same name will have to be set under Job->CustomAsimov in order to use it) and "sample" is the sample this GHOST sample will supersede |
| SeparateGammas               | if set to TRUE, the sample will not contribute to the overall gamma factors for MC stat, but a separate set of them will be added for this sample (through the SHAPE systematic technology); NB: you need to re-run at least the "b" step if you want to decorrelate the gammas on existing inputs (wf is not enough) |
| CorrelateGammasInRegions     | to be used only together with SeparateGammas; can be used to correlate MC stat across regions; example: "SR1:SR2,CR1:CR2:CR3" will use the same NP for the MC stat in each bin of SR1 and SR2 and in each bin of CR1, CR2 and CR3 |
| CorrelateGammasWithSample    | to be used only together with SeparateGammas; can be used to correlate MC stat of this sample with those of another sample (example usecase: when one sample is derived from another one through reweighting, and they are both used in the fit) |
| Morphing                     | add this to each template you have, to do a template fit / morphing; syntax is `<name-of-parameter>,<value-corresponding-to-this-template>`; the POI should be set to `<name-of-parameter>` |
| BuildPullTable               | if set to TRUE or NORM-ONLY, create tables showing the post-fit acceptance effect of nuisance parameter pulls for this sample, set to NORM+SHAPE to include the bin-by-bin effect |
| MCstatScale                  | scales up/down the MC stat size; useful to project sensitivity to larger MC samples |
| SystFromSample               | set it to TRUE (default FALSE) to have this sample inheriting all the systematics from another sample |
