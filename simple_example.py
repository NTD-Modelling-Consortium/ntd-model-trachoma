from trachoma.trachoma_functions import *
import multiprocessing
import time
from joblib import Parallel, delayed
num_cores = multiprocessing.cpu_count()


coverageFileName = 'scen2a.xlsx'

params = {'N': 2500,
          'av_I_duration' : 2,
          'av_ID_duration':200/7,
          'inf_red':0.45,
          'min_ID':11, #Parameters relating to duration of infection period, including ID period
          'av_D_duration':300/7,
          'dis_red':0.3,
          'min_D':1, #Parameters relating to duration of disease period
          'v_1':1,'v_2':2.6,
          'phi':1.4,
          'epsilon':0.5,#Parameters relating to lambda function- calculating force of infection
          'ag':0.0016, #Decay rate for age component of decline in D duration
          'prop_age':0.4, #Proportion of immunity attributed to age
          #Parameters relating to MDA
          'MDA_Cov':0.8,
          'MDA_Eff': 0.85, # Efficacy of treatment
          'rho':0.3,
          'nweeks_year':52,
          'babiesMaxAge':0.5, #Note this is years, need to check it converts to weeks later
          'youngChildMaxAge':9,#Note this is years, need to check it converts to weeks later
          'olderChildMaxAge':15, #Note this is years, need to check it converts to weeks later
          'b1':1,#this relates to bacterial load function
          'ep2':0.114,
          'n_inf_sev':30}

                       


sim_params = {'timesim':4787, 
              'burnin': 3900,
              'N_MDA':5,
              'nsim':10}


demog = {'tau': 0.0004807692, 
         'max_age': 3120,
         'mean_age': 1040}
previous_rounds = 0
vals = Set_inits(params, demog, sim_params)
vals = Seed_infection(params, vals)


MDAData = readCoverageData(coverageFileName)
TX_Mat = Tx_matrix_2(MDAData, params, previous_rounds)
MDA_dates = getMDADates(MDAData)



Start_date = date(2014, 1, 1)
End_date = date(2030,12,31)
MDA_times = get_MDA_times(MDA_dates, Start_date, sim_params['burnin'])
MDA_dates = [date(2015,5,1),date(2016,5,1),date(2017,5,1),date(2018,5,1),date(2019,5,1),date(2020,5,1),date(2021,5,1)]

MDA_times = get_MDA_times(MDA_dates, Start_date, sim_params['burnin'])
sim_params = {'timesim':5307, 
                   'burnin': 3900,
                   'N_MDA':len(MDA_times),
                   'nsim':10}


# year to make endgame output
outputYear = range(2017, 2041)
outputTimes = getOutputTimes(outputYear)
outputTimes = get_MDA_times(outputTimes, Start_date, sim_params['burnin'])
#############################################################################################################################
######################################################################################################################################################

# do sims
MDAData = readCoverageData(coverageFileName = 'scen1.xlsx')

MDA_dates = getMDADates(MDAData)
MDA_times = get_MDA_times(MDA_dates, Start_date, sim_params['burnin'])



numSims = 15
print( f'Running {numSims} simulations on {num_cores} cores' )


# randomly pick indices for number of simulations
# indices = np.random.choice(a=range(200), size = numSims, replace=False)
# indices = range(200)
start_time = time.time()

beta = 0.14
vals = Set_inits(params, demog, sim_params)
vals = Seed_infection(params, vals)
Tx_mat = Tx_matrix_2(MDAData, params, 0)
# run simulations in parallel
results = Parallel(n_jobs=num_cores)(
         delayed(sim_Ind_MDA_Include_Survey)(params=params, Tx_mat = Tx_mat, 
                                             vals = vals, timesim = sim_params['timesim'],
                                             demog=demog, bet=beta, MDA_times = MDA_times, 
                                             MDAData=MDAData,moutputTimes= outputTimes, 
                                             seed = i) for i in range(numSims))





outs = getResults(results, demog, params, outputYear)
outs
outs.to_csv('outs1.csv',index=False)

