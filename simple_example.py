import numpy as np

from trachoma.trachoma_functions import *
import multiprocessing
import time
from joblib import Parallel, delayed
num_cores = multiprocessing.cpu_count()
import pickle


#############################################################################################################################
#############################################################################################################################

# set seed and save state, leave seed=None for random data, or a value like seed=0 for consistent run-to-run data
numpy_seed = None
np.random.seed(numpy_seed)
numpy_state = np.random.get_state()

# initialize parameters, sim_params, and demography

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
          'babiesMaxAge':0.5, #Note this is years, need to check it converts to weeks later
          'youngChildMaxAge':9,#Note this is years, need to check it converts to weeks later
          'olderChildMaxAge':15, #Note this is years, need to check it converts to weeks later
          'b1':1,#this relates to bacterial load function
          'ep2':0.114,
          'n_inf_sev':38,
          'TestSensitivity': 0.96,
          'TestSpecificity': 0.965,
          'SecularTrendIndicator': 0,
          'SecularTrendYearlyBetaDecrease': 0.05,
          'vacc_prob_block_transmission':  0.8, 
          'vacc_reduce_bacterial_load': 0.5, 
          'vacc_reduce_duration': 0.5,  
          'vacc_waning_length': 52 * 5}


sim_params = {'timesim':52*23, 
              'burnin': 26,
              'N_MDA':5,
              'nsim':10}


demog = {'tau': 0.0004807692, 
         'max_age': 3120,
         'mean_age': 1040}



previous_rounds = 0


Start_date = date(2019,1, 1)
End_date = date(2030,12,31)
#############################################################################################################################
#############################################################################################################################
# import pickle file and beta values. 

# begin by specifying the name of the IU
IUID = 'BDI06375'
# path to pickle file
PickleFilePath = 'OutputVals_' + IUID + '.p'
# load pickle file
pickleData =  pickle.load(open(PickleFilePath, 'rb'))
# load beta file
BetaFilePath = 'InputBet_' + IUID + '.csv'
# read in parameters
allBetas = pd.read_csv(BetaFilePath)
#############################################################################################################################
#############################################################################################################################
# make sure the N parameter is the same as the number of people in the pickle file
a = pickleData[1]
params['N'] = len(a['IndI'])
#############################################################################################################################
#############################################################################################################################
# which years to make endgame output specify and convert these to simulation time
outputYear = range(2019, 2041)
outputTimes = getOutputTimes(outputYear)
outputTimes = get_Intervention_times(outputTimes, Start_date, sim_params['burnin'])


#############################################################################################################################
#############################################################################################################################

# generate MDA data from coverage file
scenario = '2c'
coverageFileName = 'scen' + scenario + '.csv'
MDAData = readPlatformData(coverageFileName, "MDA")
MDA_dates = getInterventionDates(MDAData)
MDA_times = get_Intervention_times(MDA_dates, Start_date, sim_params['burnin'])
sim_params['N_MDA'] = len(MDA_times)

VaccData = readPlatformData(coverageFileName, "Vaccine")
Vaccine_dates = getInterventionDates(VaccData)
vacc_times = get_Intervention_times(Vaccine_dates, Start_date, sim_params['burnin'])
sim_params['N_Vaccines'] = len(vacc_times)
#############################################################################################################################
#############################################################################################################################

# decide how many sims we will run
numSims = 2
print( f'Running {numSims} simulations on {num_cores} cores' )
start = time.time()

# generate seed
# we set the seed for generating the numpy states below, leave seed=None for random data, or a value like seed=0 for consistent run-to-run data
seed = None
np.random.seed(seed)
# we generate a numpy state for each simulation by saving a state. If the seed is set above, this will be consistent from run to run
numpy_states = list(map(lambda s: seed_to_state(s), np.random.randint(2^32, size=numSims)))

#############################################################################################################################
#############################################################################################################################
# run as many simulations as specified
results = Parallel(n_jobs=num_cores)(
         delayed(run_single_simulation)(pickleData = pickleData[i], 
                                        params = params, 
                                        timesim = sim_params['timesim'],
                                        burnin = sim_params['burnin'],
                                        demog=demog, 
                                        beta = allBetas.beta[i], 
                                        MDA_times = MDA_times, 
                                        MDAData=MDAData, 
                                        vacc_times = vacc_times, 
                                        VaccData = VaccData,
                                        outputTimes= outputTimes, 
                                        index = i,
                                        numpy_state=numpy_states[i]) for i in range(numSims))


print(time.time()- start)
#############################################################################################################################
#############################################################################################################################
# collate and output IHME data

outsIHME = getResultsIHME(results, demog, params, outputYear)
outsIHME.to_csv('outsIHME'+ scenario +'.csv',index=False)



#############################################################################################################################
#############################################################################################################################
# collate and output IPM data
MDAAgeRanges = getInterventionAgeRanges(coverageFileName, "MDA")
VaccAgeRanges = getInterventionAgeRanges(coverageFileName, "Vaccine")
outsIPM = getResultsIPM(results, demog, params, outputYear, MDAAgeRanges, VaccAgeRanges)
outsIPM.to_csv('outsIPM'+ scenario +'.csv',index=False)

