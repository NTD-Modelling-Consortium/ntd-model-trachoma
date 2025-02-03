from unittest.mock import MagicMock
import numpy as np
import copy
import unittest
import numpy.testing as npt
from ntdmc_trachoma.trachoma_functions import *

class TestMDAFunctionality(unittest.TestCase):
    # start by defining parameters for the run
    def setUp(self):
        self.params = {'N': 5000,
                       'av_I_duration': 2,
                       'av_ID_duration': 200/7,
                       'inf_red': 0.45,
                       'min_ID': 11, 
                       'av_D_duration': 300/7,
                       'min_D': 1, 
                       'dis_red':0.3,
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
                       'vacc_waning_length': 52 * 5,
                       'importation_rate': 0,
                       'importation_reduction_rate': 1,
                       'infection_risk_shape' : 1}

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
                        [2019.0, 0, 0.5, 0.8, 1, 2],
                        [2020.0, 1, 10, 0.8, 1, 2]]
        # pick some times corresponding to these MDA's. This isn't really important for the tests
        # which use this set of MDA_times as they are considered separately, we test the 
        # MDA works correctly when done concurrently on different age groups with different coverage
        # using MDA_times_concurrent and MDA_times_concurrent
        self.MDA_times = np.array([5200, 5252, 5304])

        self.MDADataConcurrent = [[2018.0, 0, 100.0, 0.1, 0, 2],
                    [2018.0, 0, 0.5, 0.4, 1, 2]]
        # adjust MDA's so that they occur at the same time to test that this is accounted for correctly
        self.MDA_times_concurrent = np.array([5200, 5200])

        seed = 0
        np.random.seed(seed)
        numpy_states = list(map(lambda s: self.seed_to_state(s), np.random.randint(2**32, size=1)))
        
        # we initiate a set of values for the starting point and will also make it so everyone is infected so that
        # we can easily test the MDA's due to everyone who is cured having a change in IndI and bact_load variables which we can check against
        self.vals = Set_inits(params=self.params, demog=self.demog, sim_params=self.sim_params, MDAData=self.MDAData, numpy_state=numpy_states[0])  
        self.vals['IndI'] = np.ones(self.params['N'])
        self.vals['No_Inf'] = np.ones(self.params['N'])
        self.vals['bact_load'] = bacterialLoad(params=self.params, vals=self.vals)
        # set a number of repetitions for the MDA so that we are likely to reach the level of tolerance we want
        # as with only one MDA there will be enough randomness to possibly be outside of this tolerance range
        self.nReps = 100

    def seed_to_state(self, seed):
        np.random.seed(seed)
        return np.random.get_state()

    def testMDAOnBabies(self):
        MDA_round = np.where(self.MDA_times == self.MDA_times[1])[0]
        # set up propCured, which is the proportion of people who are cured each MDA. This is later used
        # to calculate the mean proportion of people cured over multiple repetitions
        propCured = []
        for _ in range(self.nReps):
            valsTest = copy.deepcopy(self.vals)
            # set everyone's age to 0 so that they all count as babies
            valsTest['Age'] = np.zeros(self.params['N'])
            # keep track of how many people are infected before the MDA.
            preMDAInf = sum(valsTest['IndI'])
            for l in range(len(MDA_round)):
                MDA_round_current = MDA_round[l]
                # we want to get the data corresponding to this MDA from the MDAdata
                ageStart, ageEnd, cov, label, systematic_non_compliance = get_MDA_params(self.MDAData, MDA_round_current, valsTest)
                # if cov or systematic non compliance have changed we need to re-draw the treatment probabilities
                # check if these have changed here, and if they have, then we re-draw the probabilities
                valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
                # do the MDA for the age range specified by ageStart and ageEnd
                valsTest, _ = MDA_timestep_Age_range(valsTest, self.params, ageStart, ageEnd, 0, label, self.demog)
                # keep track of how many people are infected after the MDA.
                # this way we can calculate the mean proportion of people cured
                postMDAInf = sum(valsTest['IndI'])
                propCured.append((preMDAInf - postMDAInf) / preMDAInf)
        # calculated the expected proportion of people cured via the efficacy of the MDA along with the coverage
        # for babies, this is multiplies by 0.5 as the drug is modelled as being less effective for babies
        expectedProportionCured = 0.5 * self.params['MDA_Eff'] * cov
        npt.assert_allclose(np.mean(propCured), expectedProportionCured , atol=5e-03, err_msg="The values are not close enough")


    def testMDAOnNonBabies(self):
        MDA_round = np.where(self.MDA_times == self.MDA_times[0])[0]
        # set up propCured, which is the proportion of people who are cured each MDA. This is later used
        # to calculate the mean proportion of people cured over multiple repetitions
        propCured = []
        for _ in range(self.nReps):
            valsTest = copy.deepcopy(self.vals)
            # set everyone's age to 20 so that they all count as non-babies
            valsTest['Age'] =  20 * np.ones(self.params['N']) * 52
            # keep track of how many people are infected before the MDA.
            preMDAInf = sum(valsTest['IndI'])
            for l in range(len(MDA_round)):
                MDA_round_current = MDA_round[l]
                # we want to get the data corresponding to this MDA from the MDAdata
                ageStart, ageEnd, cov, label, systematic_non_compliance = get_MDA_params(self.MDAData, MDA_round_current, valsTest)
                # if cov or systematic non compliance have changed we need to re-draw the treatment probabilities
                # check if these have changed here, and if they have, then we re-draw the probabilities
                valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
                # do the MDA for the age range specified by ageStart and ageEnd
                valsTest, _ = MDA_timestep_Age_range(valsTest, self.params, ageStart, ageEnd, 0, label, self.demog)
                # keep track of how many people are infected after the MDA.
                # this way we can calculate the mean proportion of people cured
                postMDAInf = sum(valsTest['IndI'])
                propCured.append((preMDAInf - postMDAInf) / preMDAInf)
        # calculated the expected proportion of people cured via the efficacy of the MDA along with the coverage
        expectedProportionCured = self.params['MDA_Eff'] * cov
        npt.assert_allclose(np.mean(propCured), expectedProportionCured, atol=5e-03, err_msg="The values are not close enough")

    # tests so far have checked that for a population all aged the same, the MDA works correctly
    # the next test ensures that when the MDA is done on a typical population, no one outside of the specified
    # range of ages given by ageStart and ageEnd is given MDA
    def testAgeRangeForMDAInTypicalPopulation(self):
        MDA_round = np.where(self.MDA_times == self.MDA_times[2])[0]
        # set up propCured, which is the proportion of people who are cured each MDA. This is later used
        # to calculate the mean proportion of people cured over multiple repetitions
        propCured = []
        for _ in range(self.nReps):
            valsTest = copy.deepcopy(self.vals)
            for l in range(len(MDA_round)):
                MDA_round_current = MDA_round[l]
                # we want to get the data corresponding to this MDA from the MDAdata
                ageStart, ageEnd, cov, label, systematic_non_compliance = get_MDA_params(self.MDAData, MDA_round_current, valsTest)
                # if cov or systematic non compliance have changed we need to re-draw the treatment probabilities
                # check if these have changed here, and if they have, then we re-draw the probabilities
                valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
                # test the MDA function directly so that we can see the people who are targeted and cured by the MDA
                curedPeople, treatedPeople = doMDAAgeRange(valsTest, self.params, ageStart, ageEnd)
                # get the ages of people who were treated and cured
                treatedAges = valsTest['Age'][treatedPeople.astype(int)]
                curedAges = valsTest['Age'][curedPeople.astype(int)]

        # check that the ages of people who are treated and cured are within the bounds of ageStart and ageEnd
        assert(max(treatedAges)/52 <= ageEnd)
        assert(min(treatedAges)/52 >= ageStart)
        assert(max(curedAges)/52 <= ageEnd)
        assert(max(curedAges)/52 >= ageStart)
        

    def testNotTreatingOutsideOfTargetAgeRange(self):
        # this test is testing that when we target only babies and set the whole population
        # to be non-babies that no one will be treated
        MDA_round = np.where(self.MDA_times == self.MDA_times[1])[0]
        # set up propCured, which is the proportion of people who are cured each MDA. This is later used
        # to calculate the mean proportion of people cured over multiple repetitions
        propCured = []
        for _ in range(self.nReps):
            valsTest = copy.deepcopy(self.vals)
            # set everyone's age to 20 so that they all count as non-babies
            valsTest['Age'] = 20 * np.ones(self.params['N']) * 52
            preMDAInf = sum(valsTest['IndI'])
            for l in range(len(MDA_round)):
                MDA_round_current = MDA_round[l]
                # we want to get the data corresponding to this MDA from the MDAdata
                ageStart, ageEnd, cov, label, systematic_non_compliance = get_MDA_params(self.MDAData, MDA_round_current, valsTest)
                # if cov or systematic non compliance have changed we need to re-draw the treatment probabilities
                # check if these have changed here, and if they have, then we re-draw the probabilities
                valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
                # do the MDA for the age range specified by ageStart and ageEnd
                valsTest, _ = MDA_timestep_Age_range(valsTest, self.params, ageStart, ageEnd, 0, label, self.demog)
                # keep track of how many people are infected after the MDA.
                # this way we can calculate the mean proportion of people cured
                postMDAInf = sum(valsTest['IndI'])
                propCured.append((preMDAInf - postMDAInf) / preMDAInf)
        # calculated the expected proportion of people cured.
        # this should be 0 in this case, as everyone in the population is outside the target age range
        expectedProportionCured = 0
        npt.assert_allclose(np.mean(propCured), expectedProportionCured, atol=0, err_msg="The values are not close enough")


        
    # this is to test that if we have multiple MDA's in the same time point then we will do them appropriately
    # the test is again for babies and non-babies with different coverages for the tests for babies and non-babies
    # this also tests that the re-drawing of the treatment probabilities is working correctly,
    # at least in terms of the mean of the probabilities after re-drawing them (it doesn't test the rank correlation of this)
    def testDoingConcurrentMDAs(self):
        MDA_round = np.where(self.MDA_times_concurrent == self.MDA_times_concurrent[0])[0]
        # we will track the proportion of babies cured and the proportion of non-babies cured
        propCuredBabies = []
        propCuredNonBabies = []
        for _ in range(self.nReps):
            for l in range(len(MDA_round)):
                valsTest = copy.deepcopy(self.vals)
                valsTest['IndI'] = np.ones(self.params['N'])
                preMDAInf = sum(valsTest['IndI'])
                if(l == 0):
                    # when we target non-babies make everyone a non-baby
                    valsTest['Age'] = 20 * np.ones(self.params['N']) * 52
                else:
                    # when we target babies make everyone a baby
                    valsTest['Age'] = 0.1 * np.ones(self.params['N']) * 52
                MDA_round_current = MDA_round[l]
                # we want to get the data corresponding to this MDA from the MDAdata
                ageStart, ageEnd, cov, label, systematic_non_compliance = get_MDA_params(self.MDADataConcurrent, MDA_round_current, valsTest)
                # if cov or systematic non compliance have changed we need to re-draw the treatment probabilities
                # check if these have changed here, and if they have, then we re-draw the probabilities
                valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
                # do the MDA for the age range specified by ageStart and ageEnd
                valsTest, _ = MDA_timestep_Age_range(valsTest, self.params, ageStart, ageEnd, 0, label, self.demog)
                # keep track of how many people are infected after the MDA.
                # this way we can calculate the mean proportion of people cured
                postMDAInf = sum(valsTest['IndI'])
                if(l == 0):
                    propCuredNonBabies.append(((preMDAInf - postMDAInf) / preMDAInf))
                else:
                    propCuredBabies.append(((preMDAInf - postMDAInf) / preMDAInf))
        # calculated the expected proportion of people cured via the efficacy of the MDA along with the coverage
        # for babies, this is multiplies by 0.5 as the drug is modelled as being less effective for babies
        expectedProportionCuredNonBabies = self.params['MDA_Eff'] * self.MDADataConcurrent[0][3]
        expectedProportionCuredBabies = self.params['MDA_Eff'] * self.MDADataConcurrent[1][3] * 0.5

        npt.assert_allclose(np.mean(propCuredNonBabies), expectedProportionCuredNonBabies, atol=5e-03, err_msg="The values are not close enough for non babies")
        npt.assert_allclose(np.mean(propCuredBabies), expectedProportionCuredBabies, atol=5e-03, err_msg="The values are not close enough for babies")


    # this test is to check that once we re-draw the treatment probabilities people are still in the same rank order
    # e.g. if you were the most likely to get treated with the previous probabilities, you still are after re-drawing them
    def testRankCorellationOfTreatmentProbabilites(self):
        valsTest = copy.deepcopy(self.vals)
        # store the rank of the treatment probabilities as we will compare this to an updated set of probabilites
        rank1 = np.argsort(valsTest['treatProbability'])
        # we also store the mean of treatment probabilities to compare with the expected mean
        mean1 = np.mean(valsTest['treatProbability'])
        expected_mean1 = valsTest["MDA_coverage"] 
        # ensure that the coverage is changed so that we will have to re-draw the probabilites
        if( valsTest["MDA_coverage"] == 0):
            cov = 0.5
        else:
            cov = valsTest["MDA_coverage"] * valsTest["MDA_coverage"]
        systematic_non_compliance = 0.3
        # re-draw the probabilites
        valsTest = check_if_we_need_to_redraw_probability_of_treatment(cov, systematic_non_compliance, valsTest)
        # store the updated rank of the treatment probabilities as we will compare this to original rank
        rank2 = np.argsort(valsTest['treatProbability'])
        # we also store the mean of treatment probabilities to compare with the expected mean
        mean2 = np.mean(valsTest['treatProbability'])
        expected_mean2 = cov

        npt.assert_equal(actual = rank1, desired = rank2)
        npt.assert_allclose(mean1, expected_mean1, atol=2e-03, err_msg="The values are not close enough")
        npt.assert_allclose(mean2, expected_mean2, atol=2e-03, err_msg="The values are not close enough")


    def test_doMDAAgeRange_on_babies_needs_random_below_half(self):
        np.random.uniform = MagicMock(side_effect=[
            np.array([0] * 2), # uniform result for treating baies
            np.array([0.6, 0.4]), # uniform result for treatment working
            np.array([]), # uniform result for treating older
            np.array([]) # uniform result for older workign
        ])
        vals = {
            'Age' : np.zeros(2), # 2 zero-week olds
            'treatProbability': np.ones(2), # garanteeed to be treated
        }
        params = { 'MDA_Eff': 1.0}
        curedPeople, treatedPeople = doMDAAgeRange(vals, params, ageStart=0, ageEnd=0.5)
        npt.assert_array_equal(curedPeople, [1])
        npt.assert_array_equal(treatedPeople, [0, 1])

    # Additional tests
    # treat probability for baby excludes baby
    # 26 week old baby is treated with half efficacy 
    # ageStart = 12/52 then babies under 12 weeks no treated (I think there is a minor bug here)
    # For all tests, need to test for both whether the AgeStart is less than 0.5 or not
    # Treat probability of adult excludes people
    # MDA_Eff effects treated older people

    # Tests for MDA_timestep_age_range
    # No one cured leaves everyone with original indi
    # Someon cured sets IndI to zero, but leaves other infected untouched
    # Someon cured sets bact_load to zero, but leaves other infected untouched