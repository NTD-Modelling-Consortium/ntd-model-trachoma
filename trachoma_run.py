from trachoma.trachoma_simulations import *

# Decide how long simulation you want and when you want MDA to be carried out
sim_params = dict(
    timesim=12 * 52,  # years total duration of simulation (*52 so in weeks) including burn-in
    burnin=0 * 52,   # years burn-in to be discarded (*52 so in weeks)
    Freq_MDA=1,       # 1 for annual, 0.5 for biannual, 2 for 2 years etc.
    N_MDA=13,         # no. of rounds of MDA to be carried out in simulation
    n_sim=100         # number of simulations
)

# General parameters relating to transmission of infection
parameters = dict(

    # Population size
    N=1000,

    # Parameters relating to duration of infection period, including ID period
    av_I_duration=2,
    av_ID_duration=200 / 7,
    inf_red=0.45,
    min_ID=11,

    # Parameters relating to duration of disease period
    av_D_duration=300 / 7,
    dis_red=0.3,
    min_D=1,

    # Parameters relating to lambda function - calculating force of infection
    v_1=1,
    v_2=2.6,
    phi=1.4,
    epsilon=0.5,

    # Parameters relating to MDA
    MDA_Cov=0.8,   # MDA coverage
    MDA_Eff=0.85,  # Efficacy of treatment
    rho=0.3        # Correlation parameter for systematic non-compliance function
)

# Demography parameters
demog = dict(
    tau=1 / (40 * 52),  # death rate in weeks^-1
    max_age=60 * 52,    # maximum age in population
    mean_age=20 * 52    # mean age in population
)

# run the simulations
results = Trachoma_Simulation(parameters=parameters, sim_params=sim_params, demog=demog)

# save the results
results.to_json('trachoma_results.json')
