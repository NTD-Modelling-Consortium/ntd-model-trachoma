import numpy as np
import numpy.testing as npt
import unittest  

from trachoma.trachoma_functions import * 

class TestGetLambdaStep(unittest.TestCase): 

    def setUp(self):
        self.bet = 1
        self.timesim = 100 * 52
        self.burnin = 50 * 52

        self.params = {
            'SecularTrendIndicator': 1,
            'SecularTrendYearlyBetaDecrease': 0,
            'nSims':1
        }

    def test_WhenSpecifyBetaAsOneValue(self):
        bet = 0.6
        betas = calculateWeeklyBetas(self.timesim, self.burnin, bet, self.params)
        
        # Assert array is as expected
        npt.assert_array_almost_equal(len(betas), self.timesim)
        npt.assert_array_almost_equal(betas[0], bet)
        npt.assert_array_almost_equal(betas[-1], bet)  

    def test_WhenSpecifyBetaAsOneValueWithSecularTrend(self):
        bet = 0.6
        self.params['SecularTrendYearlyBetaDecrease'] = 0.01
        betas = calculateWeeklyBetas(self.timesim, self.burnin, bet, self.params)
        
        # Assert array is as expected
        npt.assert_array_almost_equal(len(betas), self.timesim)
        npt.assert_array_almost_equal(betas[0], bet)
        npt.assert_array_almost_equal(betas[-1], bet * (1 - self.params['SecularTrendYearlyBetaDecrease']) ** (self.timesim/52 - self.burnin/52))  


    def test_WhenSpecifyYearlyBeta(self):
        betas = np.zeros(int(self.timesim/52))
        start_beta = 0.1
        end_beta = 0.06
        
        for i in range(len(betas)):
            betas[i] = start_beta + (end_beta - start_beta) * i/(len(betas)-1)
        
        
        betas = calculateWeeklyBetas(self.timesim, self.burnin, betas, self.params)
        
        npt.assert_array_almost_equal(len(betas), self.timesim)
        npt.assert_array_almost_equal(betas[0], start_beta)
        npt.assert_array_almost_equal(betas[-1], end_beta)  
        npt.assert_array_almost_equal(betas[52], start_beta + (end_beta - start_beta) * 1/(int(self.timesim/52)-1))  
        npt.assert_array_almost_equal(betas[520], start_beta + (end_beta - start_beta) * 10/(int(self.timesim/52)-1))  

    def test_invalid_bet_length(self):
        invalid_bet = np.ones(25)  # Invalid length, not 1, timesim/52, or timesim
        with self.assertRaises(ValueError) as context:
            calculateWeeklyBetas(self.timesim, self.burnin, invalid_bet, self.params)
        
        # Check if the error message contains expected text
        self.assertIn("Invalid length for bet", str(context.exception))