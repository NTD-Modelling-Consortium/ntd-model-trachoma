import numpy as np
import copy
import unittest
import numpy.testing as npt
from trachoma.trachoma_functions import *

class TestMDAFunctionality(unittest.TestCase):
    # start by defining parameters for the run
    def setUp(self):
        self.params = {'N': 2500,
                       'av_I_duration': 2,
                       'av_ID_duration': 200/7,
                       'inf_red': 0.45,
                       'min_ID': 11, 
                       'av_D_duration': 300/7,
                       'min_D': 1, 
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
                       'SecularTrendYearlyBetaDecrease': 0.01,
                       'vacc_prob_block_transmission': 0, 
                       'vacc_reduce_bacterial_load': 0, 
                       'vacc_reduce_duration': 0,
                       'vacc_coverage': 0,  
                       'vacc_waning_length': 52 * 5}

        burnin = 100 * 52
        timesim = burnin + 21 * 52

        self.sim_params = {'timesim': timesim, 
                           'burnin': burnin,
                           'N_MDA': 5,
                           'n_sim': 16}

        self.demog = {'tau': 0.0004807692, 
                      'max_age': 3120,
                      'mean_age': 1040}

        # we set up two MDA's, one targeting everyone and another targeting just babies
        # as there is a difference in the efficacy of the drug for babies, so we want to
        # check this is working ok
        self.MDAData = [[2018.0, 0, 100.0, 0.1, 0, 2],
                        [2019.0, 0, 0.5, 0.8, 1, 2]]
        # pick some corresponding to these MDA's. This isn't really important for this test
        self.MDA_times = np.array([5200, 5252])

        seed = None
        np.random.seed(seed)
        numpy_states = list(map(lambda s: self.seed_to_state(s), np.random.randint(2**32, size=1)))
        
        # we initiate a set of values for the starting point and will also make it so everyone is infected so that
        # we can easily test the MDA's due to everyone who is cured having a change in IndI and bact_load variables which we can check against
        self.vals = Set_inits(params=self.params, demog=self.demog, sim_params=self.sim_params, MDAData=self.MDAData, numpy_state=numpy_states[0])  
        self.vals['IndI'] = np.ones(self.params['N'])
        self.vals['No_Inf'] = np.ones(self.params['N'])
        self.vals['bact_load'] = bacterialLoad(range(self.params['N']), params=self.params, vals=self.vals)
    
    def seed_to_state(self, seed):
        np.random.seed(seed)
        return np.random.get_state()

    def testMDAOnBabies(self):
        MDA_round = np.where(self.MDA_times == self.MDA_times[1])[0]
        propCured = []
        bactLoadReduction = []
        for _ in range(100):
            valsTest = copy.deepcopy(self.vals)
            # set everyone's age to 0 so that they all count as babies
            valsTest['Age'] = np.zeros(self.params['N'])
            preMDAInf = sum(valsTest['IndI'])
            preMDAbactLoad = sum(valsTest['bact_load'])
            for l in range(len(MDA_round)):
                MDA_round_current = MDA_round[l]
                ageStart, ageEnd, cov, systematic_non_compliance = get_MDA_params(self.MDAData, MDA_round_current, valsTest)
                valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
                out = MDA_timestep_Age_range(valsTest, self.params, ageStart, ageEnd)
                valsTest = out[0]
                postMDAInf = sum(valsTest['IndI'])
                postMDAbactLoad = sum(valsTest['bact_load'])
                propCured.append((preMDAInf - postMDAInf) / preMDAInf)
                bactLoadReduction.append((preMDAbactLoad - postMDAbactLoad) / preMDAbactLoad)
        
        npt.assert_allclose(np.mean(propCured), 0.5 * self.params['MDA_Eff'] * cov, atol=5e-03, err_msg="The values are not close enough")
        npt.assert_allclose(np.mean(bactLoadReduction), 0.5 * self.params['MDA_Eff'] * cov, atol=5e-03, err_msg="The values are not close enough")

    def testMDAOnAdults(self):
        MDA_round = np.where(self.MDA_times == self.MDA_times[0])[0]
        propCured = []
        bactLoadReduction = []
        for _ in range(100):
            valsTest = copy.deepcopy(self.vals)
            valsTest['Age'] = 20 * np.ones(self.params['N']) * 52
            preMDAInf = sum(valsTest['IndI'])
            preMDAbactLoad = sum(valsTest['bact_load'])
            for l in range(len(MDA_round)):
                MDA_round_current = MDA_round[l]
                ageStart, ageEnd, cov, systematic_non_compliance = get_MDA_params(self.MDAData, MDA_round_current, valsTest)
                valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
                out = MDA_timestep_Age_range(valsTest, self.params, ageStart, ageEnd)
                valsTest = out[0]
                postMDAInf = sum(valsTest['IndI'])
                postMDAbactLoad = sum(valsTest['bact_load'])
                propCured.append((preMDAInf - postMDAInf) / preMDAInf)
                bactLoadReduction.append((preMDAbactLoad - postMDAbactLoad) / preMDAbactLoad)
        
        npt.assert_allclose(np.mean(propCured), self.params['MDA_Eff'] * cov, atol=5e-03, err_msg="The values are not close enough")
        npt.assert_allclose(np.mean(bactLoadReduction), self.params['MDA_Eff'] * cov, atol=5e-03, err_msg="The values are not close enough")

    def testNotTreatingOutsideOfTargetAgeRange(self):
        MDA_round = np.where(self.MDA_times == self.MDA_times[1])[0]
        propCured = []
        for _ in range(100):
            valsTest = copy.deepcopy(self.vals)
            # set everyone's age to 20 so that none count as babies
            valsTest['Age'] = 20 * np.ones(self.params['N']) * 52
            preMDAInf = sum(valsTest['IndI'])
            for l in range(len(MDA_round)):
                MDA_round_current = MDA_round[l]
                ageStart, ageEnd, cov, systematic_non_compliance = get_MDA_params(self.MDAData, MDA_round_current, valsTest)
                valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
                out = MDA_timestep_Age_range(valsTest, self.params, ageStart, ageEnd)
                valsTest = out[0]
                postMDAInf = sum(valsTest['IndI'])
                propCured.append((preMDAInf - postMDAInf) / preMDAInf)
        
        npt.assert_allclose(np.mean(propCured), 0, atol=0, err_msg="The values are not close enough")

if __name__ == '__main__':
    unittest.main()
