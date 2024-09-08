# Trachoma Simulation Model

To run the trachoma simulation model, import the `Trachoma_Simulation()` function from the `trachoma_simulations` module in the `trachoma` package.

The `Trachoma_Simulation()` function requires the following inputs.

    BetFilePath: str
        This is the path to the input CSV file with the
        random seed and beta to be used for each simulation.

    MDAFilePath: str
        This is the path to the input CSV file with the
        first and last year of the simulations and with
        the first and last year of MDA.

    PrevFilePath: str
        This is the path where the output CSV file with
        the simulated prevalence will be saved.

    SaveOutput: bool
        If True, the last state of the simulations will
        be saved in a pickle file. If False, the last
        state of the simulations will not be saved.

    OutSimFilePath: str
        This is the path where the output pickle file with
        the last state of the simulations will be saved. It
        is only required when SaveOutput = True.

    InSimFilePath: str
        This is the path where the input pickle file with
        the last state of the simulations has been saved.
        If this is provided, the code will skip the burnin
        and resume the previous simulations from this state.
        If this is not provided, the code will start new
        simulations from scratch, including the burnin.

Numerous different examples are included in the Jupyter notebook `trachoma_tests.ipynb` in the `tests` folder.

### How to run

Install `pipenv` and `setuptools` according to the instructions for your OS, e.g. 
```commandline
pip install pipenv
pip install setuptools
```

`cd` to the project directory and run:
```commandline
pipenv install .
```

At this point you will need to add files to the root of the project directory. Please ask your 
administrator for these files.

Then you can run an example file using `pipenv`:
```commandline
pipenv run python simple_example.py
```

or open up the `pipenv` shell and run:
```commandline
pipenv shell
python simple_example.py
```

### Notes about Python versions

If you have both python 2 and 3 installed, you may need to provide `pipenv` commands with the correct version,
e.g. `pipenv install . --python 3`, `pipenv run python3 simple_example.py`, and `pipenv shell; python3 simple_example.py`
along with using `pip3`.

Note that the project has deprecated dependencies that require Python 3.8, so if you have a newer version
installed and are using an IDE, you may need to set the Python interpreter to 3.8.

### How to test

There is currently one working test. It can be run using the command

```commandline
cd tests
python -m unittest test_endtoend.py
```

### Building the docs

You'll need to have [the sphinx static site generator](https://www.sphinx-doc.org) installed.  A good way to install Sphinx is to use [`pipx`](pipx.pypa.io).

From the root of the repository:

```shell
sphinx-build -b html docs/source docs/build/html
```

You can now visualise the documentation webiste by opening `docs/build/html/index.html` with your web browser.