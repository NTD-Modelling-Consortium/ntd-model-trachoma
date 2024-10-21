import numpy as np
import numpy.testing as npt
import unittest  # Import unittest

from trachoma.trachoma_functions import * 

class TestGetLambdaStep(unittest.TestCase): 

    def setUp(self):
        # Setup common parameters
        self.demog = {'tau': 0.0004807692, 
                    'max_age': 3120,
                    'mean_age': 1040}
        self.bet = 1
        self.params = {'N': 2,
                    'ep2':0.114,
                    'vacc_prob_block_transmission': 0.5,
                    'vacc_reduce_bacterial_load': 0.5, 
                    'vacc_waning_length': 52 * 5}

        self.vals = dict(
            IndD=np.zeros(self.params['N']), 
            vaccinated=np.full(self.params['N'], fill_value=False, dtype=bool),
            time_since_vaccinated=np.zeros(self.params['N']),
            No_Inf = np.ones(self.params['N']), 
            T_ID = 10*np.ones(self.params['N']), 
        )
        
    
    def test_WhenAllHaveActiveInfection(self):
        # First Bacterial load calculation
        bact_load = bacterialLoad(self.params, self.vals)
        # Assert array is as expected
        npt.assert_array_almost_equal(np.array([1,1]),
                                    bact_load)
    
    def test_WhenNotActivelyInfected(self):
        # Modify vals and recalculate bacterialLoad
        self.vals['T_ID'][0] = 0
        
        bact_load = bacterialLoad(self.params, self.vals)
        # Assert array is as expected
        npt.assert_array_almost_equal(np.array([0,1]),
                                    bact_load)
        
    def test_When2ndInfection(self):
        # Change vaccinated and time_since_vaccinated and recalculate lambda_step
        self.vals['T_ID'][0] = 10
        self.vals['No_Inf'][0] = 2

        bact_load = bacterialLoad(self.params, self.vals)
        # Assert array is as expected
        npt.assert_array_almost_equal(np.array([1*np.exp(-self.params['ep2']),1]),
                                    bact_load)
        
    def test_WhenInfectedButVaccinated(self):
        # Change 
        self.vals['vaccinated'][0] = True

        bact_load = bacterialLoad(self.params, self.vals)
        # Assert array is as expected
        npt.assert_array_almost_equal(np.array([1 * self.params["vacc_reduce_bacterial_load"],1]),
                                    bact_load)
        
    def test_When10thInfectionButVaccinated(self):
        # Change 
        numberOfInfection = 10
        self.vals['No_Inf'][0] = numberOfInfection
        self.vals['vaccinated'][0] = True

        bact_load = bacterialLoad(self.params, self.vals)
        # Assert array is as expected
        npt.assert_array_almost_equal(np.array([1 * np.exp(-(numberOfInfection-1)*self.params['ep2']) *
                                       self.params["vacc_reduce_bacterial_load"],1]),
                                    bact_load)
        
        