import csv
import unittest
from datetime import date

import pandas as pd
import numpy as np
import copy
import trachoma.trachoma_functions as tf
import multiprocessing
import time
from joblib import Parallel, delayed
num_cores = multiprocessing.cpu_count()
import pickle

class Check_Function_Including_Survey_Matches_Function_Without_Survey(unittest.TestCase):
    """
    This is an a test that the `sim_Ind_MDA` and `sim_Ind_MDA_Include_Survey` functions
    will give the same output as each other when provided with the same seed,
    and won't do so with different seeds. We only check the prevalence of infections and disease which is calculated. 

    We input 2 seeds, one for each function.
    When the simulations use the same seed, they will run both functions on the same starting 
    population and we expect the results to match. 
    When they use different seeds, the starting population will still be the same
    but the results are not expected to match each other.

    """

    def test_functions_match_with_same_seed(self):
        """
        Test that when given the same seed, we get the same results
        """
        self.run_simulation(0,0)
        with open('check_infections_no_survey.csv', newline='') as infections_no_survey_file:
            infections_no_survey = list(csv.reader(infections_no_survey_file))
        with open('check_infections_Include_survey.csv', newline='') as infections_Include_survey_file:
            infections_Include_survey = list(csv.reader(infections_Include_survey_file))
        self.assertListEqual(infections_no_survey, infections_Include_survey)
        with open('check_disease_no_survey.csv', newline='') as disease_no_survey_file:
            disease_no_survey = list(csv.reader(disease_no_survey_file))
        with open('check_disease_Include_survey.csv', newline='') as disease_Include_survey_file:
            disease_Include_survey = list(csv.reader(disease_Include_survey_file))
        self.assertListEqual(disease_no_survey, disease_Include_survey)

    def test_functions_do_not_match_with_different_seed(self):
        """
        Test that when given different seeds, we get different results
        """
        self.run_simulation(0,1)
        with open('check_infections_no_survey.csv', newline='') as infections_no_survey_file:
            infections_no_survey = list(csv.reader(infections_no_survey_file))
        with open('check_infections_Include_survey.csv', newline='') as infections_Include_survey_file:
            infections_Include_survey = list(csv.reader(infections_Include_survey_file))
        self.assertNotEqual(infections_no_survey, infections_Include_survey)
        with open('check_disease_no_survey.csv', newline='') as disease_no_survey_file:
            disease_no_survey = list(csv.reader(disease_no_survey_file))
        with open('check_disease_Include_survey.csv', newline='') as disease_Include_survey_file:
            disease_Include_survey = list(csv.reader(disease_Include_survey_file))
        self.assertNotEqual(disease_no_survey, disease_Include_survey)

    @staticmethod
    def run_simulation(seed_No_Survey, seed_Include_Survey):
        """
        Runs the simulation with a particular seed. A helper function for the tests.
        """

        #############################################################################################################################
        #############################################################################################################################

        # decide how many sims we will run
        numSims = 1

        # we set the seed for generating the numpy states when using the function with no survey options
        np.random.seed(seed_No_Survey)
        # we generate a numpy state for each simulation by saving a state. If the seed is set above, this will be consistent from run to run
        numpy_states_Include_Survey = list(map(lambda s: tf.seed_to_state(s), np.random.randint(2^32, size=numSims)))
        
        # we set the seed for generating the numpy states when using the function with survey options, which we will turn off
        np.random.seed(seed_Include_Survey)
        numpy_states_No_Survey = list(map(lambda s: tf.seed_to_state(s), np.random.randint(2^32, size=numSims)))
        # initialize parameters, sim_params, and demography

        np.random.seed(10)
        numpy_states_generate_population = list(map(lambda s: tf.seed_to_state(s), np.random.randint(2^32, size=numSims)))

        params = {'N': 500,
                  'av_I_duration' : 2,
                  'av_ID_duration':200/7,
                  'inf_red':0.45,
                  'min_ID':11, #Parameters relating to duration of infection period, including ID period
                  'av_D_duration':300/7,
                  'min_D':1, #Parameters relating to duration of disease period
                  'dis_red':0.3,
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
                  'vacc_waning_length': 52 * 5,
                  'importation_rate': 0,
                  'importation_reduction_rate': 1}


        sim_params = {'timesim':52*23,
                      'burnin': 26,
                      'N_MDA':5,
                      'nsim':10}


        demog = {'tau': 0.0004807692,
                 'max_age': 3120,
                 'mean_age': 1040}



        previous_rounds = 0


        Start_date = date(2019,1, 1)
        End_date = date(2029,12,31)
        #############################################################################################################################
        #############################################################################################################################
        #set a beta value

        beta = 0.2

        #############################################################################################################################
        #############################################################################################################################
        # which years to make endgame output specify and convert these to simulation time
        outputYear = range(Start_date.year, End_date.year)
        outputTimes = tf.getOutputTimes(outputYear)
        outputTimes = tf.get_Intervention_times(outputTimes, Start_date, sim_params['burnin'])


        #############################################################################################################################
        #############################################################################################################################

        # generate MDA data from coverage file
        # this is currently using something outside of the test folder which is bad
        # but readPlatformData is hardcoded for this
        scenario = '2c'
        coverageFileName = 'scen' + scenario + '.csv'
        MDAData = tf.readPlatformData(coverageFileName, "MDA")
        MDA_dates = tf.getInterventionDates(MDAData)
        MDA_times = tf.get_Intervention_times(MDA_dates, Start_date, sim_params['burnin'])
        sim_params['N_MDA'] = len(MDA_times)
        # this is currently using something outside of the test folder which is bad
        # but readPlatformData is hardcoded for this
        VaccData = tf.readPlatformData(coverageFileName, "Vaccine")
        Vaccine_dates = tf.getInterventionDates(VaccData)
        vacc_times = tf.get_Intervention_times(Vaccine_dates, Start_date, sim_params['burnin'])
        sim_params['N_Vaccines'] = len(vacc_times)
        #############################################################################################################################
        #############################################################################################################################

        print( f'Running {numSims} simulations on {num_cores} cores' )
        start = time.time()
        # seed the population with infection so that there are events occurring in the simulation
        # this means that when running with different seeds, there are differences in the simulation
        # which isn't guaranteed when there are no infections
        vals = tf.Set_inits(params=params, demog=demog, sim_params = sim_params, MDAData=MDAData, numpy_state =numpy_states_generate_population[0])    # Set initial conditions
        vals = tf.Seed_infection(params=params, vals=vals) # Seed infection
        v_No_Survey = copy.deepcopy(vals)
        v_Include_Survey = copy.deepcopy(vals)
        results_from_No_Survey_func = tf.sim_Ind_MDA(params=params,
                                                    vals = v_No_Survey, timesim = sim_params['timesim'],
                                                    burnin=sim_params['burnin'],
                                                    demog=demog, bet=beta,  MDA_times = MDA_times, 
                                                    MDAData=MDAData, vacc_times = vacc_times, 
                                                    numpy_state=numpy_states_Include_Survey[0])




        results_from_Include_Survey_func, _ = tf.sim_Ind_MDA_Include_Survey(params=params,
                                                    vals = v_Include_Survey, timesim = sim_params['timesim'],
                                                    burnin=sim_params['burnin'],
                                                    demog=demog, bet=beta,  MDA_times = MDA_times, 
                                                    MDAData=MDAData, vacc_times = vacc_times, VaccData = VaccData,
                                                    outputTimes= outputTimes, doSurvey=False, doIHMEOutput=False,
                                                    numpy_state=numpy_states_No_Survey[0])

        pd.DataFrame(results_from_No_Survey_func['True_Prev_Infection_children_1_9']).to_csv('check_infections_no_survey.csv', index = False) 
        pd.DataFrame(results_from_Include_Survey_func['True_Infections_Disease_children_1_9']).to_csv('check_infections_Include_survey.csv', index = False) 
        pd.DataFrame(results_from_No_Survey_func['True_Prev_Disease_children_1_9']).to_csv('check_disease_no_survey.csv', index = False) 
        pd.DataFrame(results_from_Include_Survey_func['True_Prev_Disease_children_1_9']).to_csv('check_disease_Include_survey.csv', index = False) 
        
        print(time.time()- start)
