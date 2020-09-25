# Trachoma Simulation Model

To run the trachoma simulation model, import the `Trachoma_Simulation()` function from the `trachoma_simulations` module in the `trachoma` package.

The `Trachoma_Simulation()` function requires the following inputs:
- `parameters`: dictionary containing the model parameters relating to the transmission of infection;
- `sim_params`: dictionary containing the simulation parameters;  
- `demog`: dictionary containing the demography parameters.

The `Trachoma_Simulation()` function returns a data frame with the following columns:

- `Time`: simulated time (in years);
- `Mean_Disease_Children`: simulated average prevalence of disease in children aged 1 to 9;
- `Mean_Infection_Children`: simulated average prevalence of infection in children aged 1 to 9.

The output data frame can be exported in several different formats; see `trachoma_results.json` for an example of the results in JSON format.

See also `trachoma_run.py` for an example of how to use the `Trachoma_Simulation()` function.






