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