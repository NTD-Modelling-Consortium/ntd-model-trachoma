# Trachoma Simulation Model

The `trachoma` package provides a function `run_single_simulation` to
simulate the dynamics of trachoma spreading.

To run the trachoma simulation model, import the
`Trachoma_Simulation()` function from the `trachoma_simulations`
module in the `trachoma` package.

```python
from trachoma import run_single_simulation

vals, results = run_single_simulation(
	pickleData,
	params,
	timesim,
	demog,
	beta,
	MDA_times,
	MDAData,
	outputTimes,
	index
)
```

## trachoma model outputs

- `vals` :: Dictionary describing the entire final state of the population.
- `results` :: List of `Result` dataclass instances
```python
@dataclass
class Result:
    time: float
    IndI: ndarray
    IndD: ndarray
    Age:ndarray
    NoInf: ndarray
    nMDA:Optional[ndarray] = None
    nMDADoses: Optional[ndarray] = None
    nSurvey: Optional[int] = None
    surveyPass: Optional[int] = None
    elimination: Optional[int] = None
    propMDA: Optional[ndarray] = None
```

## Model input

- `pickleData` :: pickled dictionary describing the initial state of the population.
```python
import pickle
pickleData =  pickle.load(OutputVals_BDI06375.p, 'rb')
```
- `params` :: Dictionary describing parameters governing the trachoma model
```python
params = {'N': 2500,
          'av_I_duration' : 2,
          'av_ID_duration':200/7,
          'inf_red':0.45,
          'min_ID':11, #Parameters relating to duration of infection period, including ID period
          'av_D_duration':300/7,
          'min_D':1, #Parameters relating to duration of disease period
          'v_1':1,
          'v_2':2.6,
          'phi':1.4,
          'epsilon':0.5,#Parameters relating to lambda function- calculating force of infection
          #Parameters relating to MDA
          'MDA_Cov':0.8,
          'MDA_Eff': 0.85, # Efficacy of treatment
          'rho':0.3,
          'nweeks_year':52,
          'babiesMaxAge':0.5, # In years
          'youngChildMaxAge':9,# In years
          'olderChildMaxAge':15, # In years
          'b1':1,#this relates to bacterial load function
          'ep2':0.114,
          'n_inf_sev':38,
          'TestSensitivity': 0.96,
          'TestSpecificity': 0.965}
```
- `timesim` :: Number of weeks to simulate for, excluding burnin time.
- `demog` :: Dictionary describing population demographics
```python
demog = {'tau': 0.0004807692,
         'max_age': 3120,
         'mean_age': 1040}
```
- `beta` :: Value for infection parameter beta.
- `MDA_times` :: List of integers describing MDA timesteps as number
  of weeks from beginning of simulation, including burnin period.
- `MDAData` :: List of lists, each internal list describing a
  particular MDA event, for instance:

 | Date (fractional) | Min age | Max age | coverage ratio |
 |-------------------|---------|---------|----------------|
 | 2016.5            | 2       | 17      | 0.5            |

- `outputTime` :: Array of integers giving times to output results for the Endgame project. The entries correspond to simulation time, which counts the number of weeks since the beginning of the simulation.
```python
outputYear = range(2020, 2041)
outputTimes = getOutputTimes(outputYear)
outputTimes = get_MDA_times(outputTimes, Start_date, sim_params['burnin'])
```
