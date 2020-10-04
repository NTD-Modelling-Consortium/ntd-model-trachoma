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

Install [pipenv](https://drive.google.com/drive/folders/1Or6lUkymYd_p031xKGZLcnTV4GYf-oYb) according to the instructions for your OS, then `cd` to the project directory and run:

```
	$ pipenv install . # sets up per-project python environment ('env')
	$ pipenv shell # starts a per-project shell using that env
	(ntd-model-trachoma) $ python tests/trachoma_run.py # runs the test model scenarios
```





