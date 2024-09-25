# NTDMC trachoma model

## Installation

We recommend that you install the model within a dedicated Python
virtual environment. See [venv — Creation of virtual environments](
https://docs.python.org/3/library/venv.html) for more information
about working with virtual environments.

The model can be installed with the standard `pip` package management
utility, specifying the URL of the project's GitHub repository
followed by the string `@v<x.y.z>` where `<x.y.z>` is the version
number.  The latest release's version number can be found [here](https://github.com/NTD-Modelling-Consortium/ntd-model-trachoma/tags).

```shell
   pip install git+https://github.com/NTD-Modelling-Consortium/ntd-model-trachoma.git@v1.0.1
```

## Using the model

Currently, the model is expected to be used through the top-level function `ntdmc_trachoma.run_single_simulation`:

```python
from ntd_trachoma import run_single_simulation

vals, results = run_single_simulation(
    params, vals, timesim, burnin, demog, beta, MDA_times, MDAData,
    vacc_times, VaccData, outputTimes, doSurvey, doIHMEOutput,
    numpy_state,
)
```

See [the docs](link) for more information about the meaning of the
different parameters and the objects returned by the function.

## Contributing

Start by cloning the repository:

```shell
git clone git@github.com:NTD-Modelling-Consortium/ntd-model-trachoma.git
```

Then install the package in editable mode.  It is recommend that you
do so within a Python virtual environment dedicated to this project.

```shell
cd ntd-model-trachoma && pip install --editable .[dev]
```

Specifying the `[dev]` suffix will trigger the installation of both
`ruff` and `pytest`.

Contributions to the Python code are expected to pass all checks
applied by the `ruff check` tool.  Before you commit changes, make
sure this is the case by running the `ruff check trachoma/` command
from the root of the repository.

Automated tests can be run with `pytest`:

```shell
cd tests && python -m pytest .
```

Note the current directory _must_ be set to the `tests` folder in
order for all the tests to pass.
