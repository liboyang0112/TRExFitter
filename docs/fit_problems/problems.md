# Frequent problems with the fit and how to treat them

Following sections describe problems that are frequently observed when running fits
Brief descriptions why the problems occur and what can be done about it are presented.
Last part provides a hands on tutorial where you can fix the problems in a toy configurations.


## Fit does not converge

Everyone who has ever ran some numerical fit knows that the fit can fail.
In MIGRAD, there are two very general reasons why the fit can fail:

- The convergence for the minimum fails - this means that the fit did not find the proper global minimum
- The estimation of the Hessian fails - even when you find a the correct minimum if the Hessian is inaccurate the fit will fail as this is needed to actually check the converge (see [Minimisation procedure](technical_aspects.md#minimisation-procedure)).

You should note that the fit can fail in different stages!
It can fail during the MIGRAD step (when trying to find the minimum) or even in MINOS step, when the minimum is found but the improved error estimation fails.

### How does it manifest

In TRExFitter, the quality of the fit is actually estimated twice during the `f` step.
Internally MIGRAD does one check before it print the fitted values to a terminal, then after the procedure we check the quality of the fit again.
If the fit fails TRExFitter prints warnings telling you why it failed (either problem with the convergence, problem with the Hessian or problem with the EDM being too high, > 0.001).
If this happens, the fit is retried with higher precision (setting in Minuit) up to two more times until it converges.
If the fit fails three times, TRExFitter prints an error telling you that the fit did not converge and prints the reason (status codes).

If the fit fail during the MINOS stage, errors at the end of the procedure will be printed directly from Minuit, telling you that lower/upper estimate of the error failed.
Be aware that in this case the code will continue, however, it will use the Hessian values as the error estimate of the failed parameter (Hessian estimate does exist as in order to run MINOS the fit had to converge).

### Why does this happen

There are multiple reasons why the fit fails to converge.
The likelihood may be ill-defined leading to infinite values and complex values, this is generally caused by wrong inputs/misconfiguration.
However, the fit can even fail for mathematically well-defined problems due to machine precision and numerical instabilities.
The MIGRAD algorithms relies heavily on the estimate of the gradients and Hessians, it works pretty well when these estimates are precise but can fail if these are inaccurate.
Functions that are non-smooth (discontinuous first derivative) or even continuous highly non-linear functions can lead to inaccurate estimate of the derivatives that could lead to the fit failure.

### What can be done about it

Several things can be tried to help the fit converge:

- Set `DebugLevel: 3` and search the log for Minuit/MIGRAD errors indicating that a `NaN` has been encountered - this generally means that there is a problem with your setup.
    - Cross-check the input samples, look into `Systematics/` plots and see if there are some problematic distributions
    - Try to remove Gammas with `MCstatThreshold: NONE`, the problems may come from MC stat uncertainty - this happens usually when you have very small contributions with large weights - Note this is not something you can use for a publication, it should only help you understand the behaviour!
- If nothing seems to be wrong with the input histograms, try to simplify the fit
    - Try to merge some small backgrounds together - this reduces the complexity of the problem significantly
    - Try to reduce the number of bins - less bins = simpler fit setup
    - Try to increase pruning thresholds - when it comes to stability, very small uncertainties are problematic as they lead to almost flat likelihood in some regions
    - All of these issues should help you identify the problem, but cannot be used blindly as a proper setup depends on physics motivation!
- Does your fit fail for Asimov?
  We actually can help with that, see [Fit does not even converge for Asimov](#fit-does-not-even-converge-for-asimov).


## Fit does not even converge for Asimov

This may sound like a surprise but fit can even fail for Asimov fits when we know where the minimum of the likelihood is and we even use it as a starting point so there is no need to converge.

### How does it manifest

The fit fails during the `f` step when fitting Asimov dataset.

### Why does this happen

In this case we can actually understand it!
The problem does not come from the convergence to the minimum (it starts in the minimum), but it fails because the Hessian matrix cannot be estimated properly (see [Minimisation procedure](technical_aspects.md#minimisation-procedure)).
The problem is that the fit does not require any steps to converge but the Hessian matrix starts from assumption that it is a unity matrix assuming uncorrelated NPs.
However, if the NPs _are correlated_ the fit cannot modify the Hessian as no iteration happens, thus the fit fails.

### What can be done about it

There is an easy solution to this problem, let the fit iterate - do not start from the minimum.
This can be done but setting the initial values for the NP to a random value via `SetRandomInitialNPval: 0.1`.
This will allow the fit to iterate and modify the Hessian matrix.

### Hands on

An example config file which has problems converging for Asimov dataset can be found in `/eos/user/t/tdado/TRExFitterTutorial/Configs/AsimovFails.config`.
Note: If you do not want to run directly from eos, copy the input folder from `/eos/user/t/tdado/TRExFitterTutorial/TTbarXsec_AsimovFailure` to your preferred location and change the `HistoPath: "/eos/user/t/tdado/TRExFitterTutorial/TTbarXsec_AsimovFailure/"` to the folder where you copied the files.
Copy the config to you TRExFitter folder and try to run it with

```bash
trex-fitter h config/AsimovFails.config
trex-fitter wf config/AsimovFails.config
```

The fit should fail during the fit stage.
Read the warnings/errors printer to terminal.
The fit cannot converge for Asimov.
Try to increase the !DebugLevel by setting `DebugLevel: 2` in the config file and repeat the `f` step.
Let us try to fix the failure by using random NP starting point.
Add these two lines to the `Fit` block

```
SetRandomInitialNPval: 0.1
SetRandomInitialNPvalSeed: 1234567
```

And rerun the fit again with `trex-fitter f config/AsimovFails.config`.
Did it help?

Now you can try to set different `SetRandomInitialNPvalSeed` to a different value.
Does this change the output? How big are the changes?
One interesting observation: Try to change the number of CPU used for the fit, change

```
NumCPU: 4
```

to

```
NumCPU: 1
```

And remove the random initial NPs values.
What happens?

Running multiple threads can decrease the precision of the result as numerical error can accumulate, but it should not have a large impact.


## Large constraints of NPs

First of all, we are talking about NPs here, NFs or gammas do not have the pre-fit uncertainty defined.
Sometimes the post-fit uncertainties of some NPs are smaller than the pre-fit uncertainties, this is called "constraint".
Observing constraints does not necessarily mean that there is something wrong with the fit but you should always be careful when you see them.
What the constraint represents is that the observed data do not allow such a large variation as was used as an input uncertainty.
This is a dangerous claim as we usually have dedicated measurements for the uncertainty estimation.
However, it is possible that you measurement is more sensitive than the dedicated measurement (e.g. JER has no flavour information but you are using ttbar dilepton final states dominated by b-jets.
The JER for b-jets may be different than flavour-inclusive JER.).
What you need to be careful about is that a constraint may be a spurious constraint, e.g. you have two regions where you have the same systematic uncertainty affecting both, in one region the uncertainty is constrained (let's assume it can be justified), but since both of the regions share the NP, it will be even constrained in the second region.
It may be possible that the second region on its own would not constrain the NP.
Thus all constraints needs to be understood,

### How does it manifest

The post-fit uncertainty (represented by the black line) is shorted than one on either of the sides.
You can check the `NuisParam` plots (so called "pull plots" or "Brazilian plots") or check the text file in `Fits/`.
Generally, when you see that _all of your NPs_ are constrained this probably means that the fit is unstable.

### Why does this happen

There may be a good physics arguments why this happens - large statistics and overconservative uncertainty.
However, it can also come from fit instabilities or not ideal setup physics-wise.

### What can be done about it

First of all, check the effect of the problematic NP in `Systematics/` folder to see if everything looks as expected.
Many times these problems come just from wrong/buggy inputs.

If the problems come from the instabilities it is possible that the post-fit values for NPs are not estimated precisely using the Hessian matrix you can

- Try to use MINOS, `UseMinos` for the problematic NPs and see if it improves anything

Use smoothing when appropriate (large fluctuations due to MC stat)

- If the systematic variations are fluctuating they may result in large shape effects and thus can be constrained because of only few bins.
  This can be avoided by smoothing the systematic variations.
- Even relatively small backgrounds could be responsible for constraints when the systematic variation has large fluctuations, check the plots in `Systematics/` folder for the problematic NPs
- More information about effects of smoothing is available [here](https://indico.cern.ch/event/761804/contributions/3160985/attachments/1733339/2802398/Defranchis_template_constraints.pdf)

If the problems are not technical, but come rather from the setup there are different things you can try to better understand the result:

- Try to decorrelate the problematic NP between two regions - this will tell which region is responsible for the constraints
- Try to decorrelate shape effect and normalisation effect for the problematic NP - this will tell you if the shape or the normalisation component is responsible for the constraint.
  *Note:* You can get pure shape component by using `DropNorm` setting in your config.
- If the large constrain comes from a shape effect, you can try to reduce number of bins which could reduce the shape component
- Some of the constraints can be reduces by introducing more NPs as correlations between the NPs generally reduce the constraints (e.g. btag_0 may be constrained as it has large normalisation component but adding Luminosity uncertainty could decrease the constraint as it may be highly correlated with btag_0)

### Hands on

Two example configs with different problems can be found in `/eos/user/t/tdado/TRExFitterTutorial/Configs/ManyConstraints.config` and `/eos/user/t/tdado/TRExFitterTutorial/Configs/Smoothing.config`.
Note: If you do not run directly from the eos, copy the input folder from `/eos/user/t/tdado/TRExFitterTutorial/TTbarXsec_Problematic/` to your preferred location and change the `HistoPath: "/eos/user/t/tdado/TRExFitterTutorial/TTbarXsec_Problematic/"` to the folder where you copied the files.

The first config file demonstrates a general problem with the fit.

Copy the config to you TRExFitter folder and try to run it with

```bash
trex-fitter h config/ManyConstraints.config
trex-fitter wf config/ManyConstraints.config
```

The terminal output from the fit looks reasonable (some systematics are constrained but most are not).
But when you look at the pull plots almost every systematic is constrained significantly.
Also the `Fits/` folder shows that this is not a plotting issue.
The difference between the terminal output and the final results comes from the fact that the Hessian matrix is calculated one more time after the fit converged to its final position and thus it can be different that then values printed on the screen.
If the values are very different this generally indicates fit instabilities.

Seeing that many constraints is worrisome.
First you should check all input histograms in the `Systematics/` folder.
Now try to run MINOS for all NPs via `UseMinos: all` and rerun the `f` step again (this will take ~2-3 minutes).

Do you see an improvement in the post-fit constraints?

Nevertheless, you should be very careful with this setup as it does not behave very stable.
Try to look into individual regions (the setup contains 2 different "regions" one aiming at semileptonically decaying top quark and the other one aim at hadronically decaying top quark).

Run the code with

```bash
# run with the first region
trex-fitter wf config/ManyConstraints.config Regions=SR_lep_elmu

# then run with the second region
trex-fitter wf config/ManyConstraints.config Regions=CR_Wmass_elmu
```

Does the fit behave differently when only one of the regions is used?

The second config file illustrates a surprising constraint of one NP.
Run the code with

```bash
trex-fitter h config/Smoothing.config
trex-fitter wf config/Smoothing.config
```

First of all, notice that the MINOS calculation fails for the estimate of the lower uncertainty.
This is already a sign that something is not right.
Now focus on NP "ISR_ljets" it has a very large constraint!

Now look at the `Systematics/` folder and look at the variations - they do not look large in terms of magnitude but they have funny looking shape.
Oh no, we forgot to use Smoothing for this variation and it is generated from an independent MC sample so it can have statistical fluctuations.

Try applying smoothing to this systematic.
Open the config file and search for "ISR_ljets" you will see it is defined twice in the config, once for each region.
Add `Smoothing: 400` to both of the blocks (the 400 tells the code that you want to smooth variations that are statistically uncorrelated - they originate from different MC samples - this applies only when the `TTBARRESONANCE` smoothing option is used).
Close the file and rerun the fit again starting from the histogram creation (you can change the name of the Job option to create a new folder)

```bash
trex-fitter h config/Smoothing.config
trex-fitter wf config/Smoothing.config
```

Now check the pulls for this NP again, did it help? Also, check other NPs, did this affect also another NPs?

You will still see many NPs that are constrained but when you look in `Systematics/` folder these uncertainties are large and thus it makes sense that these are constrained.
Nevertheless, try to run only with one region in the fit to see where the constrain comes from:

```bash
# run with the first region
trex-fitter wf config/Smoothing.config Regions=SR_lep_elmu

# then run with the second region
trex-fitter wf config/Smoothing.config Regions=CR_Wmass_elmu
```

Now remove the Smoothing for this NP again and try to run MINOS on this NP via `UseMINOS: ISR_ljets` does this help in this case?

When you have some time left, try to play with the smoothing options, change the `SmoothingOption: TTBARRESONANCE` to different options: `MAXVARIATION` (default), `COMMONTOOLSMOOTHMONOTONIC`, `COMMONTOOLSMOOTHPARABOLIC`, `KERNELRATIOUNIFORM`, `KERNELDELTAGAUSS` or `KERNELRATIOGAUSS`.
You need to rerun the `h` step always.
Look at the `Systematics/` folder, try to focus on some large uncertainty and see how it changes.
You can also run the with (`wf`)
step and see how the smoothing option affects the expected uncertainty.


## Large and/or unexpected pulls of NPs

Again, we will be talking about NPs, NFs do not have Gaussian constraint term so even something like 2 +-0.1 in the post-fit may very well be a reasonable result!

As mentioned for the constraints, having large pulls does not necessarily mean that there is something wrong with the fit.
First of all, what should be considered as a "large pull" cannot be easily defined, but generally pulls that are larger than one sigma should be considered large.
Generally amount and size of NP pulls will be very different in an analysis with 100 events and in an analysis with 100k events.
Analyses that are totally stat dominated by statistical uncertainty generally have only very few pulls (stat fluctuations can compensate relatively large differences) while for analyses that are systematically dominated more pulls are expected as statistical fluctuations allow for only very small relative changes in each bin.

So why is it dangerous to have large pulls? The problem is if the large pulls are not obvious (not compensating for a discrepancy in between data and prediction) it may result in random pulls which is an unwanted behaviour.
Sometimes pulls that have large impact (large uncertainties) are pulled even when they are not obviously compensating for discrepancies between data and prediction.

!!! tip
    Use `SystDataPlots: TRUE` to generate plots that will show you the effect of each systematic on the data/prediction ratio - this will help you to understand if systematic uncertainty compensates for a mismodelling.

### How does it manifest

Large pulls in !NuisParameters plot or in `Fits/` folder.
Note that gammas (NPs for MC statistical uncertainties) and NFs are expected to be centred around 1, while NPs are expected to be centred around 0!

### Why does this happen

How can NP be pulled if it does not compensate for discrepancy between data and prediction?

- It is possible that your samples have large statistical fluctuations and e.g. pulling the normalisation of that sample is jsut compensating for a fluctuation in few bins.
  Check `Systematics/` folder and look for the problematic NPs.
- If the pulls are not large (less than 1 sigma) but are unexpected (by looking at plots from `SystDataPlots: TRUE` in `Systematics/` folder) - this can happen if the uncertainty is large and small pulls (in terms of sigma) can lead to large effects (in terms of bin contents).
  If the effect is large it is easy for the fit to pull the NP as it can compensate for large effects (bin-yield wise) by changing the likelihood value only slightly.
  In other words, for large systematics small pull can result in large effect on the bin content, thus it can be pulled easily.

### What can be done about it

- Again, first cross-check the problematic systematics in `Systematics/` folder
- If the problem comes from large stat fluctuations:
    - You can try to merge some backgrounds together to reduce the fluctuations
    - You can try to reduce the number of bins
- If you see unexpected pulls:
    - One option is to introduce a control region that is sensitive to the problematic NP
    - You can try to play with the binning to reduce _or increase_ the shape effect that can lead to different pulls
- A general rule with pulls is: "less bins = less pulls" (but also smaller sensitivity, usually)

### Hands on

An example config file has been prepared which is very similar to a config used in an actual measurement: `/eos/user/t/tdado/TRExFitterTutorial/Configs/DileptonXsec_emu.config`.
There are no obvious problems when running the config file

```bash
trex-fitter h config/DileptonXsec_emu.config
trex-fitter wf config/DileptonXsec_emu.config
```

This is a fit to data so some pulls are expected.
Anything that stands out?
Now try to create prefit and postfit plots (this is how the problem was indeed discovered):

```bash
trex-fitter dp config/DileptonXsec_emu.config
```

Now, compare the prefit and postfit plots, something unexpected is seen?
The !SingleTop contribution changed quite significantly, that is also seen in the pulls plots for `SingleTopXsec` but nothing crazy in the pull.
Now look at the `Systematics/` folder and find `SingleTopXsec`, now you see where the problem is?
Try to fix it in the config file.

This was a trivial copy-paste/wrong config mistake but these mistakes are very frequent!
Whenever something does not look right you need to cross-check the inputs.


## Post-fit constraint larger than 1 sigma

The post-fit constraints are expected to be ideally at 1 sigma (no constraint) or in some special cases can be smaller than 1 sigma (your data allows only smaller uncertainty than the input uncertainty)
However, there is no clear interpretation for constraint _larger_ than 1 sigma.

### How does it manifest

The post-fit uncertainty of the NP is larger than 1.
You can see it on pull plot or directly in the text file in `Fits/`

### Why does this happen

Generally, there are two reason why something like this can happen.

- Fit instability - especially when you see this for multiple NPs it usually indicates problems with the convergence that are not caught by Minuit
- Fit is stable but the errors on the NP are asymmetric - Hessian by default produces uncertainties that are symmetric, an if e.g. you true uncertainty is +1.1/-0.9 the result of the Hessian estimate of the uncertainty can be 0.9, 1.1 or something in between

### What can be done about it

- You should know by know that you should first look into the `Systematics/` folder if the inputs look ok :)
- Try to use MINOS with `UseMinos` for the problematic NP, it may correct the uncertainties
- Look at the LH scan for the problematic NP via `doLHscan` - this will tell you lot about the problem, if the LH scan shows asymmetric behaviour this is probably reason for the constraint being larger than 1
- It is possible that the LH scan will show you non-smooth behaviour - then it means you fit is unstable and you need to go back to the [Fit does not converge](#fit-does-not-converge) part

### Hands on

An example config file which has one NP constraint larger than 1 sigma can be found in `/eos/user/t/tdado/TRExFitterTutorial/Configs/ConstraintLargerThanOne.config`.
Note: If you do not want to run directly from the eos, copy the input folder from `/eos/user/t/tdado/TRExFitterTutorial/TTbarXsec_Problematic/` to your preferred location and change the `HistoPath: "/eos/user/t/tdado/TRExFitterTutorial/TTbarXsec_Problematic/"` to the folder where you copied the files.
Copy the config to you TRExFitter folder and try to run it with

```bash
trex-fitter h config/ConstraintLargerThanOne.config
trex-fitter wf config/ConstraintLargerThanOne.config
```

There were some warnings at the beginning of the fit but the fit converges well.
Now look at the `bTagSF_B_0` pull in the pull plot.
The constraint is much larger than one!
Cross-check this in the `Fits/` folder to make sure this is not a plotting issue.

Indeed, this seems like a problem.
Look into `Systematics/` folder and check if there is something wrong with the inputs for this particular uncertainty.
When everything seems fine try to run MINOS on this NP to see if this is not a problem of precision.
Set `UseMinos: ttbarXsec,bTagSF_B_2` and rerun the `f` step again.

Does the constraint value change?

To make sure run also LH scan for this particular NP via `doLHscan: bTagSF_B_2` in the Fit block of your config and rerun the `f` step again.

You can see in the `Fits/` folder that the uncertainty on the NP is very slightly asymmetric but this probably didn't cause the problems with large constraints.
The problems probably come from not perfect estimation of the Hessian matrix, resulting in the poor estimation of the post-fit uncertainty on this NP.


## Non-smooth Likelihood scan curve

The likelihood scan provides very useful cross-check for many results of the fit: the fitted value and its uncertainty, obtained significance etc..
It is recommended to run LH scan at least for the POI.
How does the LH scan work? The procedure repeats the fit multiple times (30 by default) for the specified parameter (can be POI, NF, NP or gamma) by fixing the parameter to a predefined value and then running the minimisation procedure.
The obtained likelihood value (NLL) is plotted on a vertical axis.

Recently, the LH scan fitting procedure has been adjusted slightly, previously for each of the scanned points, the fit started from the same initial values.
Since commit `e62615188d8d3cdf7ac5226edae1d21959a17f2e` (17.04.2019) the initial values are kept from the previously fitted point (second point starts from values where first point converged, etc.).
The new version leads to smoother LH curves (and also reduces running time slightly).
Nevertheless, the scan can still result in a non-smooth curve.

!!! tip
    You can control the minimum, maximum and the number of scanned points via: `LHscanMin`, `LHscanMax` and `LHscanSteps`

### How does it manifest

Kinks, non-smooth curve for some parameters in `LHoodPlots/` folder.
The expectation is that the curve follow quadratic formula (this may not be true for all fits, but the curve should at least be smooth).

### Why does this happen

Assuming that the minimum was found correctly, it represents the fact that the area around the minimum does not behave quadratically or even non smooth (first derivatives are not continuous).
This can be caused by real non-smooth properties of the likelihood (like using absolute values) or because it is numerically unstable (very large slopes numerically leading to infinities in derivatives).

### What can be done about it

First of all, consider if this is a real problem for you, if the instabilities occur in an area that is e.g. 10 sigma away from the minimum, that is probably not a problem.
If you see instabilities in regions close to the minimum, this is something to be worried about.
There ar e several thing that can cause this behaviour and, unfortunately, you need to test few setups:

- Simplify the fit by looking into individual regions, if possible
- Simplify the fit by merging some of the backgrounds together
- Simplify the fit by increasing pruning thresholds
- Simplify the fit by reducing the number of bins

All of these instructions are just to identify the problem, they should not be used for the final result!


## Unexpected high (anti-)correlations of NPs

The correlations play a vital role in the whole process of the PL fit and are responsible for most of the reduction in the post-fit uncertainty.
Thus, having correct correlations is very important!
The correlation matrix is automatically plotted for every fit in `f` step.
For every fit, check the correlation matrix and focus on large (anti-)correlations.
Some of the correlations can be easily explained (e.g. signal strength can be highly correlated with luminosity uncertainty).
Focus on the systematics with large (anti-)correlations that you did not expect.

!!! tip
    Use `CorrelationThreshold: XX` to plot only systematics with large (anti-)correlations.

### How does it manifest

Large (anti-)correlations in the correlation matrix (cold/hot spots on the heat map) for uncertainties where the (anti-)correlations are not expected

### Why does this happen

As mentioned in [Minimisation procedure](technical_aspects.md#minimisation-procedure), the correlations of NPs are obtained directly from the fitting procedure, from the Hessian matrix.
If the Hessian matrix is inaccurate, the estimation of correlation will be inaccurate.
Note that not all of the correlations that are unexpected are wrong, they may be surprising at first, but may be real (in sense that they do not come from some fit instability).

### What can be done about it

Check `Systematics/` folder for the systematics with unexpected (anti-) correlations.
Look at the plots and see if the shapes of the systematic variations look very similar.
If the variations look similar, this means the correlation is real.
If it is still problematic to have the large (anti-)correlations, you can try to increase the number of bins to further separate the shape effect of the problematic NPs.

If the variation with high (anti-) correlation look very different, this is a sign of a problem with the fit.
Unfortunately, there is no "MINOS equivalent" for Hessian, so you cannot get a more precise estimate of the correlations from the fit.
You can try to simplify the fit to understand the correlations, similarly as in [Non-smooth LIkelihood scan curve](#non-smooth-likelihood-scan-curve).
