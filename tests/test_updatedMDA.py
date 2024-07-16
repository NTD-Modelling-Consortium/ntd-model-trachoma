from trachoma.trachoma_functions import *
import numpy as np
import copy
import unittest
import numpy.testing as npt

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
          'rho':0,
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
          'SecularTrendYearlyBetaDecrease': 0.01,
          'vacc_prob_block_transmission':  0, 
          'vacc_reduce_bacterial_load': 0, 
          'vacc_reduce_duration': 0,
          'vacc_coverage': 0,  
          'vacc_waning_length': 52 * 5}

def seed_to_state(seed):
    np.random.seed(seed)
    return np.random.get_state()

burnin = 100*52
timesim = burnin + 21*52

sim_params = {'timesim': timesim, 
              'burnin': burnin,
              'N_MDA':5,
              'n_sim':16}


demog = {'tau': 0.0004807692, 
         'max_age': 3120,
         'mean_age': 1040}


# we want to test that MDA is working as expected on babies and all other people
# set up one MDA with coverage 10% for all people and one for only babies with coverage 80%
MDAData = [[2018.0, 0, 100.0, 0.1, 0, 2],
            [2019.0, 0, 0.5, 0.8, 1, 2]]
# we also need to specify the times for these MDA's so that later we can access them correctly
MDA_times = np.array([5200, 5252])

seed = None
np.random.seed(seed)
# we generate a numpy state for each simulation by saving a state. If the seed is set above, this will be consistent from run to run
numpy_states = list(map(lambda s: seed_to_state(s), np.random.randint(2^32, size=1)))
vals = Set_inits(params=params, demog=demog, sim_params = sim_params, MDAData=MDAData, numpy_state=numpy_states[0])    # Set initial conditions
# we start by assigning everyone as infected
vals['IndI'] = np.ones(params['N'])
vals['No_Inf'] = np.ones(params['N'])
vals['bact_load'] = bacterialLoad(range(params['N']),params = params, vals = vals)

def testMDAOnChildren(vals, params, MDA_times, MDAData, nMDAs, i):
    MDA_round = np.where(MDA_times == i)[0]
    propCured = []
    bactLoadReduction = []
    for _ in range(nMDAs):
        valsTest = copy.deepcopy(vals)
        valsTest['Age'] = np.zeros(params['N'])
        preMDAInf = sum(valsTest['IndI'])
        preMDAbactLoad = sum(valsTest['bact_load'])
        for l in range(len(MDA_round)):
            MDA_round_current = MDA_round[l]
            ageStart, ageEnd, cov, systematic_non_compliance = get_MDA_params(MDAData, MDA_round_current, valsTest)
            valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
            out = MDA_timestep_Age_range(valsTest, params, ageStart, ageEnd)
            valsTest = out[0]
            postMDAInf = sum(valsTest['IndI'])
            postMDAbactLoad = sum(valsTest['bact_load'])
            propCured.append((preMDAInf - postMDAInf)/preMDAInf)
            bactLoadReduction.append((preMDAbactLoad - postMDAbactLoad)/preMDAbactLoad)
    npt.assert_allclose(np.mean(propCured), 0.5*params['MDA_Eff'] * cov, atol=5e-03, err_msg="The values are not close enough")
    npt.assert_allclose(np.mean(bactLoadReduction), 0.5*params['MDA_Eff'] * cov, atol=5e-03, err_msg="The values are not close enough")




def testMDAOnAdults(vals, params, MDA_times, MDAData, nMDAs, i):    
    MDA_round = np.where(MDA_times == i)[0]
    propCured = []
    bactLoadReduction = []
    for _ in range(nMDAs):
        valsTest = copy.deepcopy(vals)
        valsTest['Age'] = 20 * np.ones(params['N']) * 52
        preMDAInf = sum(valsTest['IndI'])
        preMDAbactLoad = sum(valsTest['bact_load'])
        for l in range(len(MDA_round)):
            MDA_round_current = MDA_round[l]
            ageStart, ageEnd, cov, systematic_non_compliance = get_MDA_params(MDAData, MDA_round_current, valsTest)
            valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
            out = MDA_timestep_Age_range(valsTest, params, ageStart, ageEnd)
            valsTest = out[0]
            postMDAInf = sum(valsTest['IndI'])
            postMDAbactLoad = sum(valsTest['bact_load'])
            propCured.append((preMDAInf - postMDAInf)/preMDAInf)
            bactLoadReduction.append((preMDAbactLoad - postMDAbactLoad)/preMDAbactLoad)
    npt.assert_allclose(np.mean(propCured), params['MDA_Eff'] * cov, atol=5e-03, err_msg="The values are not close enough")
    npt.assert_allclose(np.mean(bactLoadReduction), params['MDA_Eff'] * cov, atol=5e-03, err_msg="The values are not close enough")

# this test is to ensure we don't treat outside of the target age range. 
# this is done by assigning everyone's age to be 20 years old and then 
# ensuring that no one is treated when we only target children
def testNotTreatingOutsideOfTargetAgeRange(vals, params, MDA_times, MDAData, nMDAs, i):
    MDA_round = np.where(MDA_times == i)[0]
    propCured = []
    for _ in range(nMDAs):
        valsTest = copy.deepcopy(vals)
        valsTest['Age'] = 20 * np.ones(params['N']) * 52
        preMDAInf = sum(valsTest['IndI'])
        for l in range(len(MDA_round)):
            MDA_round_current = MDA_round[l]
            ageStart, ageEnd, cov, systematic_non_compliance = get_MDA_params(MDAData, MDA_round_current, valsTest)
            valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
            out = MDA_timestep_Age_range(valsTest, params, ageStart, ageEnd)
            valsTest = out[0]
            postMDAInf = sum(valsTest['IndI'])
            propCured.append((preMDAInf - postMDAInf)/preMDAInf)
    npt.assert_allclose(np.mean(propCured), 0, atol=0, err_msg="The values are not close enough")



nMDAs = 100

testMDAOnAdults(vals, params, MDA_times, MDAData, nMDAs, MDA_times[0])

testMDAOnChildren(vals, params, MDA_times, MDAData, nMDAs, MDA_times[1])
testNotTreatingOutsideOfTargetAgeRange(vals, params, MDA_times, MDAData, nMDAs, MDA_times[1])
