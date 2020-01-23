# TRExFitter   [![build status](https://gitlab.cern.ch/TRExStats/TRExFitter/badges/master/build.svg "build status")](https://gitlab.cern.ch/TRExStats/TRExFitter/commits/master)

This package provides a framework to perform profile likelihood fits. In addition to that, many convenient features are available. TRExFitter was previously also known as TtHFitter. Here are a few important references to make use of:

* [TRExFitter twiki page](https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TtHFitter) for additional documentation and many references to further details
* [TRExFitter JIRA](https://its.cern.ch/jira/projects/TTHFITTER/summary>) (sign up to the mailing list in case you cannot access the JIRA)
* TRExFitter mailing list: [atlas-phys-stat-tthfitter](https://e-groups.cern.ch/e-groups/EgroupsSubscription.do?egroupName=atlas-phys-stat-tthfitter)
* Make sure to read the [FAQ](#faq) section and take a look at the [FitProblems Tutorial twiki](https://twiki.cern.ch/twiki/bin/view/AtlasProtected/FitProblemsTutorial) which describes common issues

Contributions to TRExFitter are welcome.
Please have a look at [CONTRIBUTING.md](CONTRIBUTING.md) to get started.



## Table of Contents
1.  [Getting the code](#getting-the-code)
2.  [Setup](#setup)
    * [Build it yourself](#build-it-yourself)
        * [Setup the code inside ATLAS environments](#setup-the-code-inside-atlas-environments)
    * [Using a Docker image](#using-a-docker-image)
        * [Setup using Docker image with Singularity](#setup-using-docker-image-with-singularity)
        * [Setup using Docker image with Docker](#setup-using-docker-image-with-docker)
3.  [How to](#how-to)
4.  [Config File](#config-file)
    * [Available settings](#available-settings)
5.  [Command line options](#command-line-options)
6.  [Ranking Plot](#ranking-plot)
7.  [Grouped Impact](#grouped-impact)
8.  [Multi-Fit](#multi-fit)
    * [Available multi-fit settings](#available-multi-fit-settings)
9.  [Running Unfolding](#running-unfolding)
10. [Input File Merging with hupdate](#input-file-merging-with-hupdate)
11. [Output Directories Structure](#output-directories-structure)
12. [ShapeFactor example](#shapefactor-example)
13. [Replacement file](#replacement-file)
14. [FAQ](#faq)
15. [TRExFitter package authors](#trexfitter-package-authors)



## Getting the code
To get the code, use the following command:
```
git clone ssh://git@gitlab.cern.ch:7999/TRExStats/TRExFitter.git
```
To get a specific tag, do the following:
```
cd TRExFitter && git checkout <tag number> && cd -
```
If you just want to run TRExFitter, you can also use a Docker image. See section [Using a docker image](#using-a-docker-image).



## Setup
There are multiple ways to setup TRExFitter. You can compile the code by yourself, or use the provided Docker images.

### Build it yourself
To setup just use the script (from any location):
```
  source setup.sh
```
(should work on any machine with access to cvmfs - provided that nothing else is set-up previously)
Note that this will setup centOS7 version of ROOT and you need centOS7 compatible OS!
If you want to setup the slc6 version of ROOT do:
```
  source setup.sh slc6
```


To compile:

1. Create a new build directory with `mkdir build && cd build`
2. Run cmake with `cmake ../`
3. Compile the code with `cmake --build ./`
4. The binary file will appear in `bin/` directory

(this will take as main code the file `util/trex-fitter.cc`)

The setup script also adds a path to the binary into your PATH and you can execute the code with `trex-fitter`

**IMPORTANT!** For the first time use you need to type `git submodule init` followed by `git submodule update`.
Every time the submodules change, you need to run `git submodule update`.

**Tip:** To recompile the code, directly from the main directory, a simple command is:
```
    cd build/ && cmake --build ./ && cd ../
```
or simply using the alias defined in the `setup.sh` script:
```
    trex-make
```

#### Setup the code inside ATLAS environments
To setup the code inside ATLAS environments such as `AnalysisBase` or `AnalysisTop`, do:
```
mkdir source build
cd source
asetup AnalysisBase,21.2.70,here
git clone --recursive ssh://git@gitlab.cern.ch:7999/TRExStats/TRExFitter.git
cd ../build/
cmake ../source && make
source  x86_64-centos7-gcc62-opt/setup.sh
```
Then, the setup is ready to execute `trex-fitter`.

### Using a Docker image
TRExFitter can be run via the provided Docker images. The image tagged `latest` corresponds to the current code version in the master branch. There are also tagged versions of the images corresponding to the TRExFitter tags, starting after version `TtHFitter-00-04-06`. An overview of the available images can be found in the [container registry](https://gitlab.cern.ch/TRExStats/TRExFitter/container_registry).

In order to use the images, you need to get a gitlab token. The token is created in the gitlab user settings, under Access Tokens ([direct link](https://gitlab.cern.ch/profile/personal_access_tokens)). Enter a name for the token (such as `docker_token`, and set the scope to `read_registry`. The token will then be shown after you clicked create.

#### Setup using Docker image with Singularity
Follow these steps to use the image via Singularity, for example on lxplus. Export the token and your username into environment variables:
```
export SINGULARITY_DOCKER_USERNAME=<CERN-username>
export SINGULARITY_DOCKER_PASSWORD=<gitlab-token>
```
Now you can run the following command:
```
singularity run --contain -B /tmp --pwd ${PWD} docker://gitlab-registry.cern.ch/trexstats/trexfitter:latest
```
Replace `latest` by another tag to get the corresponding version of the code. In the container you will directly have the `trex-fitter` executable. If you cannot see your local folders, you might need to mount them via the `-B` flag. The TRExFitter code is located in the folder `/TRExFitter/source/TRExFitter` within the container.

You might want to add the flag ```--silent``` to suppress the warnings occuring for AnalysisBase images based on CentOS7.

#### Setup using Docker image with Docker
These steps describe how to use the image with the Docker software, for example on your own local machine. Get started by creating an account and downloading docker here: [https://hub.docker.com/](https://hub.docker.com/).

To get access to the images, start with an authentication by running the following command:
```
docker login gitlab-registry.cern.ch -u <CERN-username> -p <gitlab-token>
```
This uses the token created as described above. The docker container can then be obtained and run with the following commands:
```
docker pull gitlab-registry.cern.ch/trexstats/trexfitter:latest
docker run -it gitlab-registry.cern.ch/trexstats/trexfitter:latest
```
Replace `latest` by a tag to get a specific version of the code.



## How to
To run the code, after compiling (see [Setup](#setup)), use the command:
```
trex-fitter <action(s)> <config file> [<options>]
```
The configuration file (`<config file>`) is a text file containing all the information on the definition of samples and fit regions, including all the fit and draw options.
By default, the file  `config/myFit.config`  is loaded.
See the section [Config File](#config-file) for more details.
Take a look at `config/myFit.config` or `config/ttH2015.config` to see some example config files.
Most of the time, the only file the user has to modify to obtain their desired results is the configuration file.

The only mandatory argument, `<action(s)>`, tells TRExFitter which operation(s) to perform.
The possible operations are defined in the main file (e.g. `util/trex-fitter.C`).
For instance, if you use the default file `util/trex-fitter.C`, the available actions are:

| **Option** | **Action** |
| ---------- | ---------- |
| `u` | read efficiencies, migration/response matrices an acceptances for unfolding and then fold them |
| `h` | read input histograms (valid only if the proper option is specified in the config file) |
| `n` | read input ntuples (valid only if the proper option is specified in the config file) |
| `w` | create the RooStats xmls and workspace |
| `f` | fit the workspace |
| `l` | calculate exclusion limit |
| `s` | calculate significance |
| `d` | draw pre-fit plots |
| `p` | draw post-fit plots |
| `a` | draw separation plots |
| `r` | draw ranking plot (see [Ranking Plot](#ranking-plot)) |
| `b` | re-run smoothing (in the future also rebinning) |
| `m` | multi-fit (see [Multi-Fit](#multi-fit)) |
| `i` | grouped impact evaluation (see [Grouped Impact](#grouped-impact)) |
| `x` | run likelihood scan only, will not produce the standard fit output like pulls/correlation matrix/etc (useful with "LHscan" command line option for parallelization) |

New optional argument: `<options>`.
It is a string (so make sure to use " or ' to enclose the string if you use more than one option) defining a list of options, in the form:
```
"<option1>=<value1>,<value2>,...:<option2>=..."
```
See the section [Command line options](#command-line-options) below.



## Config File
The structure of the file should be the following:
```
<ObjectType>: <ObjectName>
  <ObjectProperty>: <Value>
  <ObjectProperty>: <Value>
  ...

<ObjectType>: <ObjectName>
  <ObjectProperty>: <Value>
  <ObjectProperty>: <Value>
  ...

...
```
NB: note the **blank** line between the objects!

The file should contain:
  * exactly one object of type `Job`
  * exactly one object of type `Fit`
  * exactly one object of type `Limit`
  * at least one object of type `Sample`
  * at least one object of type `Region`
  * any number of objects of type `Systematic` (even 0 is ok)
  * any number of objects of type `NormFactor` (even 0 is ok)

In case of unfolding you need:
  * exactly one object of type `Job`
  * exactly one object of type `Unfolding`
  * at least one object of type `TruthSample`
  * at least one object of type `UnfoldingSample`

Note that each object should have unique `<ObjectName>`.

At the beginning of TRExFitter execution, the config file used will be checked against a reference file.
The reference files for single and multi-fits are `jobSchema.config` and `multiFitSchema.config`, respectively.
These files specify which settings are allowed per block, and how the arguments should look like.
The available blocks are:
- `Job`
- `Fit`
- `Limit`
- `Significance`
- `Options`
- `Region`
- `Sample`
- `NormFactor`
- `ShapeFactor`
- `Systematic`
- `TruthSample`
- `Unfolding`
- `UnfoldingSample`
- `UnfoldingSystematic`

### Available settings
For each object type (or "block"), you can find the available settings in [our documentation (docs/Settings.md)](docs/Settings.md#standard-fit).



## Command line options
Currently the supported options are:

| **Option** | **Effect** |
| ---------- | ---------- |
| **Job**               | to to provide a new name for the output folder |
| **Regions**           | to limit the regions to use to the list specified |
| **Samples**           | to limit the samples to use to the list specified |
| **Systematics**       | to limit the systematics to use to the list specified |
| **Signal(s)**         | in case more than one SIGNAL sample is specified in your config file, you can specify which one you want to run on (for plots, workspace creation and fits/limits/significance); GHOST samples can be promoted to SIGNAL samples in this way, and multiple samples can be specified using `Samples:smp1,smp2` |
| **Exclude**           | to exclude certain Regions / Samples / Systematics |
| **Suffix**            | used for: plots, workspace, fit results, etc |
| **SaveSuffix**        | used for: saving histograms with a suffix (to be merged / renamed later, see [Input File Merging with hupdate](#input-file-merging-with-hupdate) section |
| **Update**            | if TRUE, the output .root file is updated, otherwise is overwrote |
| **StatOnlyFit**       | if TRUE, the same as Fit->StatOnlyFit |
| **StatOnly**          | if TRUE, no systematics nor MC stat uncertainties will be considered (equivalent to set StatOnly: TRUE in the Job block of the config), use `Systematics=NONE` instead to keep MC stat uncertainties |
| **Ranking**           | see [Ranking Plot](#ranking-plot) section |
| **FitResults**        | the specified fit results file will be used, for instance for post-fit plots (instead of the file `jobName/Fits/jobName.txt`) |
| **FitType**           | can be set to SPLUSB or BONLY to replace the option in the config file |
| **LumiScale**         | as the options in config file |
| **BootstrapIdx**      | see description of Bootstrap option in config (under Job) |
| **BootstrapSYst**     | see description of BootstrapSyst option in config (under Job) |
| **GroupedImpact**     | see [Grouped Impact](#grouped-impact) section |
| **OutputDir**         | see [Job settings](docs/Settings.md#job-block-settings) section |
| **LimitParamValue**   | see [Limit settings](docs/Settings.md#limit-block-settings) section (ParamValue) |
| **LHscan**            | set a NP/POI for the likelihood scan can be used for parallelization of the code |
| **Parallel2Dscan**    | run only slice of LH2D scan in x-direction can be used for parallelization of the code |
| **Parallel2Dscan**    | define which step of the parallelized 2D scan should be performed (has to be an integer between 0 and LHscanSteps-1) |
| **FitBlind**          | see [Fit settings](docs/Settings.md#fit-block-settings) section |
| **BlindedParameters** | see [Fit settings](docs/Settings.md#fit-block-settings) section |

Note: the wild-card `*` is supported, but only as last character.
Example:
```
trex-fitter n config/ttH2015.config 'Regions=HThad_ge6jge4b:Exclude=BTag_*'
```



## Ranking Plot
* The ranking plot can be created in one go, with just the command line argument `r` (after having run the nominal fit `f`).
* Since this can take too much time (and memory), for complicated fits it's better to run it in several steps:
   by specifying the command-line option `Ranking=<name/index>`, one can produce the txt input for the ranking only for a specific line of the ranking, i.e. for a single NP (specified either through its name or index). Once all the needed txt files are created (e.g. in parallel through batch jobs) with the option `Ranking=plot` they are merged to create the final plot.

* Examples:
```
# this runs the ranking in one go
trex-fitter r <config>
#these commands will first create the inputs for the ranking one by one and then merge them in the plot
trex-fitter r <config> Ranking=Lumi
trex-fitter r <config> Ranking=JES1
trex-fitter r <config> Ranking=ttXsec
trex-fitter r <config> Ranking=plot
```



## Grouped Impact
* The command line argument `i` is used to evaluate the combined impact of groups of nuisance parameters on the POI.
* Specify groups using the `SubCategory` option (for Systematics and NormFactors).
* Two groups are defined by default: "Gammas" (MC stat. impact) and "FullSyst" (full systematics impact with statistical component subtracted).
* The impact is calculated by performing a fit where the nuisance parameters in the group are fixed to their best-fit values, and then subtracting the resulting uncertainty on the POI in quadrature from the uncertainty from the nominal fit.
* The command line parameter `GroupedImpact` can be used to parallelize the impact calculations. If it is not specified, all existing groups are evaluated sequentially.
* The results are saved in `Fits/GroupedImpact*`.

Examples:

```
# evaluate impact of all groups sequentially
trex-fitter i <config>

# evaluate only the impact of Gammas
trex-fitter i <config> GroupedImpact="Gammas"
```

If the calculations are parallelized, combine the results by running the following at the end:

```
trex-fitter i <config> GroupedImpact="combine"
```



## Multi-Fit
The Multi-Fit functionality can be used to compare fit results or even to combine fit inputs from different configuration files / Jobs.

To use it you need a dedicated config file, with a structure similar to the usual ones. Example:
```
MultiFit: "myTopWS_multifit"
  Label: "My Label"
  Combine: FALSE
  Compare: TRUE
  CmeLabel: "13 TeV"
  LumiLabel: "85 pb^{-1}"
  ComparePOI: TRUE
  ComparePulls: TRUE
  CompareLimits: TRUE
  POIName: "SigXsecOverSM"
  POIRange: -10,30
  DataName: "obsData"
  CombineChByCh: TRUE

Fit: "CR"
  ConfigFile: config/myTopWS_CR.config
  Label: "CR-only"

Fit: "SR"
  ConfigFile: config/myTopWS_SR.config
  Label: "SR"
```
This config file can be run with the command line:
```
trex-fitter m config/myTopWS_multifit.config
```
  This will compare the fit results in terms of fitted NP, fitted POI and limits from the two config files specified. Notice that the fit and limits results have to be already available (they are not produced on the fly when running his multi-fit option).

To make a real combination, one needs to use the usual command options `w`, `f` and `l` together with the flag "Combine: TRUE" in the config above. Example:
```
trex-fitter mwf config/myTopWS_multifit.config
```
This will create a combined ws starting from the individual ws for the different regions in the two config files, and fit it.

You can also run ranking for the combined fit using
```
trex-fitter mr config/myTopWS_multifit.config
```

And, same as for the single fits, you can run on the individual NPs via
```
trex-fitter mr config/myTopWS_multifit.config Ranking="XXX"
```

### Available multi-fit settings
Find all available multi-fit settings in [our documentation (docs/Settings.md)](docs/Settings.md#multi-fit).



## Running Unfolding
To run the unfolding, you need to run steps: `u`, followed by `h`, then the other steps will work as usual
Example:
```
trex-fitter u test/configs/FitExampleUnfolding.config
trex-fitter h test/configs/FitExampleUnfolding.config
trex-fitter w test/configs/FitExampleUnfolding.config
trex-fitter f test/configs/FitExampleUnfolding.config
```


## Input File Merging with hupdate
A macro `hupdate` is included, which mimics hadd functionality, but without adding histograms if they have the same name.
This is useful for running different systematics in different steps (like different batch jobs) and then merging results afterwards.
`hupdate` is compiled automatically when using cmake. To explicitly request compilation, execute the following in the build folder:
```
make hupdate.exe
```
Example usage, combined with the usage of SaveSuffix:
```
make hupdate.exe
trex-fitter n ../config/ttH2015.config Systematics=BTag_B_NP1:SaveSuffix=_BTag_B_NP1
./build/bin/myFit.exe n ../config/ttH2015.config Exclude=BTag_B_NP1:SaveSuffix=_rest
./build/bin/hupdate.exe ../ttH2015/Histograms/ttH2015_HThad_4j2b_histos.root ttH2015/Histograms/ttH2015_HThad_4j2b_histos_rest.root ttH2015/Histograms/ttH2015_HThad_4j2b_histos_BTag_B_NP1.root
./build/bin/hupdate.exe ../ttH2015/Histograms/ttH2015_HThad_5j3b_histos_NEW.root ttH2015/Histograms/ttH2015_HThad_5j3b_histos.root ttH2015/Histograms/ttH2015_HThad_5j3b_histos_BTag_B_NP1.root
./build/bin/hupdate.exe ../ttH2015/Histograms/ttH2015_HThad_ge6jge4b_histos_NEW.root ttH2015/Histograms/ttH2015_HThad_ge6jge4b_histos.root ttH2015/Histograms/ttH2015_HThad_ge6jge4b_histos_BTag_B_NP1.root
trex-fitter dwf ../config/ttH2015.config
```



## Output Directories Structure
For each TRExFit object, a directory is created, with the same name as the Fit Name.
Inside this directory, at every step, some outputs are created, following the structure described above:

| **Folder** | **Content** |
| ---------- | ----------- |
| `Plots/`              | data/MC plots, pre- and post-fit, for all the Signal, Control and Validation regions, including the summary plots |
| `Tables/`             | tables in txt and tex format |
| `RooStats/`           | workspace(s) and the xmls |
| `Fits/`               | output from fits |
| `Limits/`             | outputs from the limit-setting code |
| `Significance/`       | outputs from the significance code |
| `Systematics/`        | plots for the syst variations |
| `Toys/`               | plots and ROOT files with pseudoexperiments output |
| `Histograms/`         | root file(s) with all the inputs |
| `LHoodPlots/`         | likelihood scan with respect to the specified parameter |
| `UnfoldingHistograms/`| folded histograms produced during `u` step |



## ShapeFactor example
* The following scripts create example histograms in `exampleDataDriven` directory and execute `trex-fitter` using `config/dataDriven.config`
* The example contains a control region and signal region with two bins. The shape of one of the background samples is estimated using the ShapeFactor:
```
python makeDataDriven.py
python runDataDrivenExample.py
```
The results are in `JobDataDriven`



## Replacement file
You can define placeholders in your config file, which are replaced with values specified in an external file, which is read at the beginning of TRExFitter execution. This requires adding an additional option into your config, as part of the Job block:
```
ReplacementFile: path/to/file.txt
```
The replacement file should have the following structure:
```
# comment
XXX_placeholder: 0.1
XXX_another_placeholder: 0.2
% also a comment
```
Note that all placeholders must start with ``XXX``. In your config file, you can then refer to the placeholders like this:
```
Sample: "ttbar"
  MCweight: XXX_placeholder
```
If you would like to ensure that the replacement works correctly, set your `DebugLevel` to a minimum value of 1 and check the output of the framework.



## FAQ
__The `n` step takes a very long time to run, how do I speed this up?__\
Run multiple jobs, with each job only processing all histograms for one region. For a region called `RegionA`, this is achieved by running
```
trex-fitter n your.config Regions="RegionA"
```
Run one job per region, then proceed with `w` only when all jobs are finished. Multiple regions per job can be processed via `Regions="RegionA,RegionB"`.

It is possible to split up the histogram creation into even more jobs, splitting up systematics or samples. See section see [Input File Merging with hupdate](#input-file-merging-with-hupdate) for more details, in this case an extra step is required before proceding with `w`.

The runtime of the `n` step also depends on the amount of events in the ntuple(s) that are being read for each region. Events in the ntuple(s) that do not pass the region selection requirements further increase the execution time. A pre-selection or splitting of the input events across several files (e.g. one file per region) can decrease the execution time.

__There is some problem with the fit. What can I do?__\
Have a look at the [TRExFitter twiki](https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TtHFitter), and also the [FitProblemsTutorial twiki page](https://twiki.cern.ch/twiki/bin/view/AtlasProtected/FitProblemsTutorial).

__Why is the E in TRExFitter capitalized?__\
TRExFitter stands for "Top Related Experiment Fitter".

__How do the smoothing algorithms work?__\
Some info can be found in these [slides](https://indico.cern.ch/event/691683/contributions/2873279/attachments/1593521/2522846/PruningSmoothing.pdf). More information about `TTBARRESONANCE` is found in these [slides](https://indico.cern.ch/event/669913/contributions/2769795/attachments/1549339/2433688/ttres-fullunblind-smooth2-summary2.pdf).
Further info regarding `TTBARRESONANCE`: The Smoothing parameter in the Systematics area can be set to 40 to treat the systematic uncertainty as correlated with the nominal (e.g. when obtained via reweighting) or 400 to treat it as uncorrelated with the nominal (e.g. for two-point systematics, then the statistical uncertainties on nominal sample and systematics variation are added in quadrature and compared to the smoothing threshold).

__I am getting warnings about underconstrained nuisance parameters. Should I worry?__\
Have a look at the relevant [part of the FitProblemsTutorial twiki](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/FitProblemsTutorial#4_5_Post_fit_constraint_larger_t). Typically an underconstraint indicates an issue, though there can be cases where it is expected. See [this notebook](https://cernbox.cern.ch/index.php/s/DiPdvBlRBQPfHEy) for such an example.

__Can I just use more bins to gain sensitivity?__\
While more bings generally increase sensitivity, it is important to keep two things in mind. Limit the size of the MC statistical uncertainties ("gammas") to at most 20% per bin. Larger values can bias the signal extraction, see [these slides](https://indico.cern.ch/event/615262/contributions/2484815/) for a study and the corresponding recommendation. Furthermore, one aspect that is not considered in the fit is the statistical uncertainty in the templates that describe the ±1σ variations of a nuisance parameter. This uncertainty increases when adding more bins, and can lead to fluctuations in the templates. The nuisance parameter then does not describe physical effects, but rather fluctuations. It is very important to closely study all templates found in the `Systematics/` folder and verify that the input distributions to the fit (solid lines) look reasonable. These plots are obtained by enabling the settings `SystControlPlots` and `SystErrorBars` (both enabled by default).

__How do the automatic binning algorithms work?__\
The algorithms for `TransfoD` and `TransfoF` can be found in [these slides](https://indico.cern.ch/event/455289/contributions/1953694), `TransfoJ` is found in [slides here](https://indico.cern.ch/event/472696/contributions/1992693/). See also [this thesis](https://cds.cern.ch/record/2296985/), section 5.3.1. In practice, the `TransfoD` algorithm was found to work well, and a popular setting for the two parameters is to have both equal to the same integer. The amount of bins in the distribution is then equal to the sum, i.e. `"AutoBin","TransfoD",4,4` will create a distribution with 8 bins.


## TRExFitter package authors
Managers:

* Michele Pinamonti [michele.pinamonti@gmail.com](michele.pinamonti@gmail.com)
* Loic Valery [loic.valery@cern.ch](loic.valery@cern.ch)

Development and support team:

* Alexander Held [alexander.held@cern.ch](alexander.held@cern.ch)
* Tomas Dado [tomas.dado@cern.ch](tomas.dado@cern.ch)
