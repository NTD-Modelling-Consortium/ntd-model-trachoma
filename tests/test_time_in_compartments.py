import numpy as np
import numpy.testing as npt
import math
import unittest  # Import unittest

from trachoma.trachoma_functions import * 

class TestTimeInCompartmentsIsSameAsRModel(unittest.TestCase): 

    # check that the way the python code does it gives the same result
    def setUp(self):
        # Setup common parameters
        self.R_params = {'nun': 1/((7.7*7)/math.log(2)),
                        'nu0' : 0.0033,
                        "inf_red" : 0.4509,
                        'nun_A' : 1/((1*7)/math.log(2)),
                        'nu0_A':0.0050,
                        'inf_red_A':0.3065}
        
        self.demog = {'tau': 0.0004807692, 
                    'max_age': 3120,
                    'mean_age': 1040}
        self.bet = 1
        self.params = {'N': 2,
                       'min_ID': (1/self.R_params['nun'])/7,
                       'inf_red': self.R_params['inf_red'],
                       'dis_red': self.R_params['inf_red_A'],
                       'min_D': (1/self.R_params['nun_A'])/7,
                       'vacc_reduce_duration': 0
                    }

        self.vals = dict(
            Ind_ID_period_base=np.ones(self.params['N']) * (1/self.R_params['nu0'])/7, 
            Ind_D_period_base = np.ones(self.params['N']) * (1/self.R_params['nu0_A'])/7,
            vaccinated=np.full(self.params['N'], fill_value=False, dtype=bool),
            No_Inf=np.ones(self.params['N']),
            Age = np.ones(self.params['N'])
        )

    def test_firstInfection(self):
        IDs = ID_period_function(range(self.params['N']), self.params, self.vals)
        R_times = (self.R_params['nu0']-self.R_params['nun'])*np.exp(-self.R_params['inf_red']*(1-1))+ self.R_params['nun']
        npt.assert_array_almost_equal(IDs, np.round(1/(7*R_times)*np.ones(self.params['N'])))

    def test_tenthInfection(self):
        self.vals['No_Inf'] = self.vals['No_Inf'] * 10
        IDs = ID_period_function(range(self.params['N']), self.params, self.vals)
        R_times = (self.R_params['nu0']-self.R_params['nun'])*np.exp(-self.R_params['inf_red']*(10-1))+ self.R_params['nun']
        npt.assert_array_almost_equal(IDs, np.round(1/(7*R_times)*np.ones(self.params['N'])))
    
    def test_firstDiseased(self):
        Ds = D_period_function(self.vals['Ind_D_period_base'], self.vals['No_Inf'], self.params, self.vals['Age'])
        R_times = (self.R_params['nu0_A'] - self.R_params['nun_A'] )* np.exp(-self.R_params['inf_red_A']*(1-1))+ self.R_params['nun_A']
        npt.assert_array_almost_equal(Ds, np.round(1/(7*R_times)*np.ones(self.params['N'])))

    def test_tenthDiseased(self):
        self.vals['No_Inf'] = self.vals['No_Inf'] * 10
        Ds = D_period_function(self.vals['Ind_D_period_base'], self.vals['No_Inf'], self.params, self.vals['Age'])
        R_times = (self.R_params['nu0_A'] - self.R_params['nun_A'] )* np.exp(-self.R_params['inf_red_A']*(10-1))+ self.R_params['nun_A']
        npt.assert_array_almost_equal(Ds, np.round(1/(7*R_times)*np.ones(self.params['N'])))