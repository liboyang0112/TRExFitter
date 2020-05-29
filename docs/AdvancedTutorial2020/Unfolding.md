# Unfolding using profile likelihood

## Basic idea

Standard profile likelihood fitting can be used for unfolding.
The key is to transform the the unfolding problem into a standard problem of fitting normalisation of distributions.
This can be done easily, by _folding_ a truth distribution using a response matrix bin-by-bin.
This procedure will result in one detector level distribution per truth bin.
Normalisations of each detector level distributions can be left free floating in the fit.
Fitting of the detector level folded distributions gives the normalisation of each truth bin, which is the desired unfolded distribution.
For more details, check this [presentation](https://indico.cern.ch/event/890060/contributions/3754199/attachments/1991168/3320058/Unfolding_with_TRExFitter.pdf).

## Folding

The folding step can be done within `TRExFitter`, currently only with the histograms as an input.
The relevant truth distribution and response matrix need to be provided for the signal sample.
Note that for each systematic uncertainty, the response matrix has to be provided to generate the templates needed for systematic variations.
Alternatively, selection efficiency (probability of a given truth even to be selected per bin), migration matrix and acceptance (probability of reco event to be reconstructed for events that do not pass the truth selection) can be provided instead of the response matrix.

!!! tip "Definition of response matrix"
    The response matrix is defined as selection efficiency times migration matrix divided by acceptance. All bin contents must be $\leq 1$.

Now, let us look at how this is done in `TRExFitter`.
We will use a simple examples that is also used in our CI tests.
The folding step is the first step that needs to be run for an unfolding job.

```bash
trex-fitter u test/configs/FitExampleUnfolding.config
```

### Job block

Let us look at what is defined in the config.
The first thing that you should notice is that in the `Job` part, new path settings appeared: `MigrationPath`, `SelectionEffPath`, `ResponseMatrixPath`, `ResponseMatrixFile`, `AcceptancePath` and `AcceptanceFile`.
These options only set paths for the relevant files and histograms.
The standard `TRExFitter` logic is used here that `Path`, `File`, `Name` are combined to get the relevant histograms.
The settings in `Job` can be overwritten by the same setting in `Region` which can in turn be overwritten by settings from  `UnfoldingSample` or `Unfolding Systematics`.

### Fit block

Another important difference to the standard config file is the `FitType: UNFOLDING` setting in the `Fit` block.
This tells the code that the main result of this configuration is unfolding and not a regular fit.

### Unfolding block
This is a completely new block used only in unfolding jobs.
The first setting `MatrixOrientation` simply tells the code what your definition of horizontal and vertical axes is.

`TruthDistributionPath`, `TruthDistributionFile`, `TruthDistributionName` tell the code where to look for the truth distribution.
Note that there is only one truth distribution, even when multiple signal regions for unfolding are defined.

`NumberOfTruthBins` just tells the code how many truth bins there are, this is mainly to allow some automatic cross-checks in the code.

`Tau` parameter allows to add a Gaussian (Tikhonov) constraint term to a given bin. The syntax is as follows: `binIndex1:value,binIndex2:value,binIndex3:value`.

An important option is `NominalTruthSample`, which tells the code which of the defined `TruthSample` is the nominal one that is used for folding (others are used just for plotting).

`UnfoldingResultMin` and `UnfoldingResultMax` set range of the normalisation parameters for the truth bins.

All the other options are just cosmetic options for the final plot.

### TruthSample

These blocks just define the truth distribution.
Multiple distributions can be defined, these will then be plotted in the final unfolded distribution.
Note that one sample needs to be defined as nominal using the `NominalTruthSample` setting in the `Unfolding` block.

### Region

The standard definition of regions is extended to define paths for the migration matrix, selection efficiency and acceptance as is shown in the example config.
Again, `NumberOfRecoBins` is used to define the number of reco bins, mostly for the sake of automatic tests.

### UnfoldingSample

A new block is provided for the signal samples that will be used in the unfolding.
The relevant paths are defined here.
Note, multiple unfolding samples can be used at the same time, this allows to combine multiple detector level regions for unfolding.

### UnfoldingSystematic

This new block is used for defining systematic uncertainties that act on the signal and need the folding step.
Similarly to standard systematics, `MigrationNameUp` can be defined (same for selection efficiency and acceptance) which will replace the nominal migration matrix when this uncertainty is being folded.

!!! tip "Corelating uncertainties"
    You can correlate the standard `Systematic` and `UnfoldingSystematic` using the same `NuisanceParameter` name as is done in the example config with `bTag_SF_B_eigen0`.

Now with the basic parts of the config file explained, have a look at the output produced by the `u` step.
First, you will see that some plots are created in the `FitExampleUnfolding` folder.
These plots show the migration matrix as well as the response matrix.
The folded distributions can be found in `FitExampleUnfolding/UnfoldingHistograms/`.
You can inspect them.

## Fitting
Now that all the hard work is done to prepare the inputs needed for unfolding, you can run the usual steps:

```bash
trex-fitter h test/configs/FitExampleUnfolding.config
```

Which will read the input histograms.
It will also read the folded distributions from the previous step without you needing to change anything!

Now we can make the workspace and run the fit

```bash
trex-fitter wf test/configs/FitExampleUnfolding.config
```

Apart from the usual plots made by `TRExFitter`, you should also see a new plot called `UnfoldedData` that shows the unfolded distributions as well as all truth distributions defined in the config.

You can also produce the standard pre/post fit plots running:

```bash
trex-fitter dp test/configs/FitExampleUnfolding.config
```

!!! question "Interpreting the results"
    Why does the postfit distribution agree almost perfectly with data even with very few uncertainties defined?
    What happens when you repeat the process, but remove the `Tau` parameter completely?

!!! tip "Running toys"
    It may be very useful to run pseudoexperiments to test the unfolding procedure. You can do this using the standard config option in `Fit` block: `FitToys: XXX` where `XXX` is the number of toys to be used. Note that for toys, the Asimov prediction will be used.
