import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing

from trachoma.trachoma_functions import *

num_cores = multiprocessing.cpu_count()

def Trachoma_Simulation(parameters, sim_params, demog):

    # Longitudinal simulations:
    # First Set initial conditions and seed infection
    vals = Set_inits(params=parameters, demog=demog)    # Set initial conditions
    vals = Seed_infection(params=parameters, vals=vals) # Seed infection

    # Set times for MDA to be carried out
    MDA_times = Set_t_MDA(sim_params=sim_params)

    # Specify which individuals get treated at each MDA, create treatment matrix
    Tx_mat = Tx_matrix(params=parameters, sim_params=sim_params)

    # Setting beta; Corresponds to roughly 20% true TF in ages 1-9. AB can provide further
    # ranges of beta values to correspond to baseline prevalence estimates
    bet = np.random.uniform(size=sim_params['n_sim'], low=0.05, high=0.12)

    # Run multiple simulations
    def multiple_simulations(i):

        np.random.seed(1 + i + 5)

        out = sim_Ind_MDA(params=parameters, vals=vals, timesim=sim_params['timesim'], demog=demog, bet=bet[i], Tx_mat=Tx_mat, MDA_times=MDA_times)

        return out

    data_store_all_sim = Parallel(n_jobs=num_cores)(delayed(multiple_simulations)(i) for i in range(sim_params['n_sim']))

    # Create empty matrices length burnin:end of simulation
    True_Prev_Infection_children_1_9 = np.zeros(shape=(sim_params['timesim'] - sim_params['burnin'], sim_params['n_sim']))
    True_Prev_Disease_children_1_9 = np.zeros(shape=(sim_params['timesim'] - sim_params['burnin'], sim_params['n_sim']))
    True_Prev_Disease = np.zeros(shape=(sim_params['timesim'] - sim_params['burnin'], sim_params['n_sim']))
    True_Prev_Infection = np.zeros(shape=(sim_params['timesim'] - sim_params['burnin'], sim_params['n_sim']))
    Time = np.arange(sim_params['timesim'] - sim_params['burnin'])

    # Fill with output
    for i in range(sim_params['n_sim']):

        True_Prev_Disease[:, i] = data_store_all_sim[i]['True_Prev_Disease'][sim_params['burnin']: sim_params['timesim']]
        True_Prev_Infection[:, i] = data_store_all_sim[i]['True_Prev_Infection'][sim_params['burnin']: sim_params['timesim']]
        True_Prev_Infection_children_1_9[:, i] = data_store_all_sim[i]['True_Prev_Infection_children_1_9'][sim_params['burnin']: sim_params['timesim']]
        True_Prev_Disease_children_1_9[:, i] = data_store_all_sim[i]['True_Prev_Disease_children_1_9'][sim_params['burnin']: sim_params['timesim']]

    # Save means for all
    results = pd.DataFrame({'Time': Time / 52,
                            'Mean_Disease_Children': np.mean(True_Prev_Disease_children_1_9, axis=1),
                            'Mean_Infection_Children': np.mean(True_Prev_Infection_children_1_9, axis=1)})

    return results

