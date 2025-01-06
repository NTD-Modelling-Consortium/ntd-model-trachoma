import numpy as np
import numpy.testing as npt
import unittest  # Import unittest

from trachoma.trachoma_functions import * 

class TestSurvey(unittest.TestCase): 

    def initialSurvey(self, surveyPrev, MDA_times, nMDAWholePop):
        numMDAForSurvey = nMDAWholePop + numMDAsBeforeNextSurvey(surveyPrev)
        surveyTime = min(MDA_times) + (5 * 52) + 25
        if surveyPrev <= 0.05:
            surveyTime = min(MDA_times) + 25
        return numMDAForSurvey, surveyTime
    
    def subsequentSurvey(self, surveyPrev, timesim, nMDAWholePop, MDA_times, i, surveyPass):
        numMDAForSurvey = 0
        surveyTime = min(MDA_times) + (5 * 52) + 25
        if surveyPrev <= 0.05:
                surveyPass += 1
        else:
            surveyPass = 0

        # if we have passed 2 surveys, we won't do another one, so set the surveytime to be after simulation ends
        if surveyPass == 2:
            surveyTime = timesim + 10
        # if we have passed 1 survey, we will do another in 2 years time
        elif surveyPass == 1:
            surveyTime = i + 104  
        else: # if we didn't pass the survey, we will do another survey after a number of MDAs based on the prevalence. 
            # Assume that these MDAs must cover a significant portion of the population so call these nMDAWholePop.
            # add the number of MDAs already done to the number of MDAs to be done before the next survey
            numMDAForSurvey = nMDAWholePop + numMDAsBeforeNextSurvey(surveyPrev)
        return numMDAForSurvey, surveyTime, surveyPass

    def create_diseased(self, initial_prevalence: float):
        vals = Set_inits(self.params, self.demog, self.sim_params, self.MDAData, np.random.get_state())
        ids = np.random.choice(
                range(self.params["N"]), int(initial_prevalence * self.params["N"]), replace=False
            )
        vals['IndD'][ids] = 1
        return vals
    
    def run_setup_for_Survey_test(self, seed, initial_prevalence=0.1):
        """
        Runs the simulation with a particular seed. A helper function for the tests.
        """
        # Number of simulations
        numSims = 1

        # Set the seed for reproducibility
        np.random.seed(seed)
        numpy_states = list(map(lambda s: seed_to_state(s), np.random.randint(2**32, size=numSims)))

        # Initialize parameters, sim_params, and demography
        self.params = {
            'N': 2500,
            'av_I_duration': 2,
            'av_ID_duration': 200 / 7,
            'inf_red': 0.45,
            'min_ID': 11,
            'av_D_duration': 300 / 7,
            'min_D': 1,
            'dis_red': 0.3,
            'v_1': 1,
            'v_2': 2.6,
            'phi': 1.4,
            'epsilon': 0.5,
            'MDA_Cov': 0.8,
            'MDA_Eff': 0.85,
            'rho': 0.3,
            'nweeks_year': 52,
            'babiesMaxAge': 0.5,
            'youngChildMaxAge': 9,
            'olderChildMaxAge': 15,
            'b1': 1,
            'ep2': 0.114,
            'n_inf_sev': 38,
            'TestSensitivity': 0.96,
            'TestSpecificity': 0.965,
            'SecularTrendIndicator': 0,
            'SecularTrendYearlyBetaDecrease': 0.05,
            'vacc_prob_block_transmission': 0.8,
            'vacc_reduce_bacterial_load': 0.5,
            'vacc_reduce_duration': 0.5,
            'vacc_waning_length': 52 * 5,
            'importation_rate': 0,
            'importation_reduction_rate': 1,
            'surveyCoverage': 0.4
        }

        self.sim_params = {
            'timesim': 52 * 23,
            'burnin': 26,
            'N_MDA': 5,
            'nsim': 10
        }

        self.demog = {
            'tau': 0.0004807692,
            'max_age': 3120,
            'mean_age': 1040
        }

        self.MDAData = [
            [2020.0, 1, 100, 0.8, 0, 2],
            [2020.5, 1, 15, 0.8, 1, 2]
        ]
        self.Start_date = date(2010,1, 1)
        self.MDA_dates = getInterventionDates(self.MDAData)
        self.MDA_times = get_Intervention_times(self.MDA_dates, self.Start_date, self.sim_params['burnin'])
        self.vals = self.create_diseased(initial_prevalence)

    def test_surveyPrev(self):
        # Set parameters for the test
        seed = 1
        initial_prevalence = 0.1

        # Run setup
        self.run_setup_for_Survey_test(seed, initial_prevalence)

        # Run the survey prevalence function and test against expected value
        prev, _ = returnSurveyPrev(
            self.vals,
            self.params['TestSensitivity'],
            self.params['TestSpecificity'],
            self.demog,
            0 / 52,
            self.params['surveyCoverage']
        )

        # Check the prevalence value
        npt.assert_array_almost_equal(prev, 0.121212)


    def test_surveyPrev_Higher_initial_prevalence(self):
        # Set parameters for the test
        seed = 2
        initial_prevalence = 0.2

        # Run setup
        self.run_setup_for_Survey_test(seed, initial_prevalence)

        # Run the survey prevalence function and test against expected value
        prev, _ = returnSurveyPrev(
            self.vals,
            self.params['TestSensitivity'],
            self.params['TestSpecificity'],
            self.demog,
            0 / 52,
            self.params['surveyCoverage']
        )

        # Check the prevalence value
        npt.assert_array_almost_equal(prev, 0.218919)

    def test_surveyPrev_Higher_initial_prevalence_higher_surveycoverage(self):
        # Set parameters for the test
        seed = 2
        initial_prevalence = 0.2

        # Run setup
        self.run_setup_for_Survey_test(seed, initial_prevalence)
        self.params['surveyCoverage'] = 0.8
        # Run the survey prevalence function and test against expected value
        prev, _ = returnSurveyPrev(
            self.vals,
            self.params['TestSensitivity'],
            self.params['TestSpecificity'],
            self.demog,
            0 / 52,
            self.params['surveyCoverage']
        )

        # Check the prevalence value
        npt.assert_array_almost_equal(prev, 0.20781)

    def test_initialSurvey(self):
        # Run setup
        self.run_setup_for_Survey_test(0, 0.1)
        surveyPrev = 0.1
        nMDAWholePop = 0
        numMDAsForSurvey, surveyTime = self.initialSurvey(surveyPrev, self.MDA_times, nMDAWholePop)
        npt.assert_array_equal(numMDAsForSurvey, nMDAWholePop + 3)
        npt.assert_array_equal(surveyTime, min(self.MDA_times) + (5 * 52) + 25)

        surveyPrev = 0.01
        nMDAWholePop = 0
        numMDAsForSurvey, surveyTime =  self.initialSurvey(surveyPrev, self.MDA_times, nMDAWholePop)
        npt.assert_array_equal(numMDAsForSurvey, nMDAWholePop + 100)
        npt.assert_array_equal(surveyTime, min(self.MDA_times) + 25)

    def test_subsequentSurvey(self):

        self.run_setup_for_Survey_test(0, 0.1)
        surveyPrev = 0.1
        nMDAWholePop = 0
        i = 0
        surveyPass = 0
        numMDAsForSurvey, surveyTime, surveyPass = self.subsequentSurvey(surveyPrev, self.sim_params['timesim'], nMDAWholePop, self.MDA_times, i, surveyPass)
        npt.assert_array_equal(numMDAsForSurvey, nMDAWholePop + 3)
        npt.assert_array_equal(surveyTime, min(self.MDA_times) + (5 * 52) + 25)
        npt.assert_array_equal(surveyPass, 0)

        surveyPrev = 0.01
        nMDAWholePop = 0
        i = 0
        surveyPass = 0
        numMDAsForSurvey, surveyTime, surveyPass = self.subsequentSurvey(surveyPrev, self.sim_params['timesim'], nMDAWholePop, self.MDA_times, i, surveyPass)
        npt.assert_array_equal(numMDAsForSurvey, 0)
        npt.assert_array_equal(surveyTime, i + 104)
        npt.assert_array_equal(surveyPass, 1)

        surveyPrev = 0.01
        nMDAWholePop = 0
        i = 0
        surveyPass = 1
        numMDAsForSurvey, surveyTime, surveyPass = self.subsequentSurvey(surveyPrev, self.sim_params['timesim'], nMDAWholePop, self.MDA_times, i, surveyPass)
        npt.assert_array_equal(numMDAsForSurvey, 0)
        npt.assert_array_equal(surveyTime, self.sim_params['timesim'] + 10)
        npt.assert_array_equal(surveyPass, 2)