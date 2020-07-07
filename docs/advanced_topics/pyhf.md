# TRExFitter to pyhf

[pyhf](https://scikit-hep.org/pyhf/) is a pyhon-based implementation of [HistFactory](https://cds.cern.ch/record/1456844/).
`HistFactory` specifies the model for binned template fits that `TRExFitter` (among many other frameworks in ATLAS) follows.
The `RooFit` workspace created during the `w` action follows the `HistFactory` standard.
You can take this workspace and hand it to a tool other than `TRExFitter` to run a fit, if that tool understands the `HistFactory` specification.

`TRExFitter` saves the workspace in two formats, one is using `.xml` files and the other is using `.root`.
`pyhf` includes a utility to read a workspace in `.xml` format and translate it to the `.json` format it uses.
The specification of this `.json` format is described in more detail [here in the documentation](https://scikit-hep.org/pyhf/likelihood.html).

## Creating a .json workspace

To create a workspace in `.json` format, we can use `pyhf xml2json PATH_TO_XML`.
It needs to be called with a path pointing to the main `.xml` file, which is of type "Combination" and includes:

```html
<!DOCTYPE Combination  SYSTEM 'HistFactorySchema.dtd'>
```

This file is located within the `RooStats` directory.
For a TRExFitter config with a `Job` block that looks like this:

```yaml
Job: "pyhf_example"
  Label: "..."
  ...
```

the file will be located at `pyhf_example/RooStats/pyhf_example.xml`.
The `RooStats` directory will contain additional `.xml` files for every region, with names `pyhf_example_REGION_NAME.xml`.
Regions are called "channels" in the `HistFactory` language.
You can inspect all those file, they form your serialized `HistFactory` model.

In this example structure, you can create a `.json` workspace by running

```bash
pyhf xml2json pyhf_example/RooStats/pyhf_example.xml > workspace.json
```

which will write the workspace into `workspace.json`.
This file can contain many thousand lines, and contains the relevant histograms for all templates directly in the file (instead of storing it in `.root` files which the `.xml` files point to).

### Running from different locations

If you want to build the workspace from a directory other than the one which includes the directory that contains all the TRExFitter output (including the`RooStats` folder), you can use the option `--basedir`.

#### Example 1

You are in a directory containing `pyhf_example`:

```bash
.
└── pyhf_example
    ├── Histograms
    │   └── pyhf_example_Signal_region_histos.root
    └── RooStats
        ├── pyhf_example.xml
        └── pyhf_example_Signal_region.xml
```

Run the following to create your workspace:

```bash
pyhf xml2json pyhf_example/RooStats/pyhf_example.xml > workspace.json
```

#### Example 2

You are somewhere else, and not in the folder containing `pyhf_example`:

```bash
.
└── directory
    └── pyhf_example
        ├── Histograms
        │   └── pyhf_example_Signal_region_histos.root
        └── RooStats
            ├── pyhf_example.xml
            └── pyhf_example_Signal_region.xml
```

Create your `.json` workspace like this:

```bash
pyhf xml2json --basedir=directory directory/pyhf_example/RooStats/pyhf_example.xml > workspace.json
```

## Running a fit with `pyhf`

Once you have the workspace in `.json` format, you can perform inference via `pyhf`.
As an example, the following python script runs an unconstrained maximum likelihood fit:

```python
import json
import pyhf

# use Minuit2 as minimizer through iminuit
pyhf.set_backend("numpy", pyhf.optimize.minuit_optimizer(verbose=True))

with open("workspace.json", "r") as f:
    spec = json.load(f)

workspace = pyhf.Workspace(spec)
model = workspace.model()
data = workspace.data(model)

# run the fit
result = pyhf.infer.mle.fit(data, model, return_uncertainties=True)

# format parameter names
def get_parameter_names(model):
    labels = []
    for parname in model.config.par_order:
        for i_par in range(model.config.param_set(parname).n_parameters):
            labels.append(
                "{}[bin_{}]".format(parname, i_par)
                if model.config.param_set(parname).n_parameters > 1
                else parname
            )
    return labels

bestfit = result[:, 0]
uncertainty = result[:, 1]
labels = get_parameter_names(model)

# show results
for i, label in enumerate(labels):
    print(f"{label}: {bestfit[i]} +/- {uncertainty[i]}")
```

You need to have [iminuit](https://iminuit.readthedocs.io/) installed for this, which is an interface to Minuit2.

You can also read a `pyhf`-compatible workspace directly from the `.xml` files, using this:

```python
from pyhf.readxml import parse

main_xml_path = "directory/pyhf_example/RooStats/pyhf_example.xml"
parent_folder = "directory/"

spec = parse(main_xml_path, parent_folder, track_progress=False)
```

## Caveats

There are a few things to keep in mind when fitting with `pyhf`, which may cause the results to differ from those obtained with `TRExFitter`:

### Constant parameters

`HistFactory` contains a built-in solution for luminosity uncertainties, which `TRExFitter` does not make use of.
The advantage of implementing luminosity uncertainties manually is the larger flexibility, which is an advantage when combining analyses.
See https://its.cern.ch/jira/browse/TTHFITTER-296 for more information on this choice.
The built-in so-called `lumi` parameter is set to constant in the main `.xml`:

```html
<ParamSetting Const="True">Lumi</ParamSetting>
```

As of `pyhf` 0.4.4, constant parameter settings are not passed on to the minimizers.
See [pull request #846](https://github.com/scikit-hep/pyhf/pull/846) and [issue #739](https://github.com/scikit-hep/pyhf/issues/739).

#### Pruning gammas

As a special case of constant parameters, the pruning of gammas for statistical uncertainties is implemented by setting parameters to constant in `HistFactory`.
It does not yet work in `pyhf`, see [issue #662](https://github.com/scikit-hep/pyhf/issues/662).

### Expressions

The `Expression` functionality for normalization factors in `TRExFitter` is currently not supported in `pyhf`.
See [issue #850](https://github.com/scikit-hep/pyhf/issues/850).

### Minimization fine-tuning

The `iminuit` package wraps Minuit2, which is also used for minimization in `TRExFitter`.
By default none of the default parameters are adjusted through `pyhf`, and for more complex fits some fine-tuning may be needed.

## Example

A minimal example of a TRExFitter setup converted to `pyhf` is shown in this repository: [alexander-held/template-fit-workflows](https://github.com/alexander-held/template-fit-workflows).
