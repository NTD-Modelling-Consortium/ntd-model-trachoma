import numpy as np
import numpy.testing as npt
import unittest  # Import unittest

from ntdmc_trachoma.trachoma_functions import *

class TestGetLambdaStep(unittest.TestCase): 

    def setUp(self):
        # Setup common parameters
        self.demog = {'tau': 0.0004807692, 
                    'max_age': 3120,
                    'mean_age': 1040}
        self.bet = 1
        self.params = {'N': 3,
                    'v_1': 1,
                    'v_2': 2.6,
                    'phi': 1.4,
                    'epsilon': 0.5,
                    'vacc_prob_block_transmission': 0.5,
                    'vacc_waning_length': 52 * 5,
                    'infection_risk_shape': 1}

        self.Age = np.array([1 * 52, 10 * 51, 50 * 52])
        self.vals = dict(
            IndD=np.zeros(self.params['N']), # no one is diseased
            bact_load=np.ones(self.params['N']),  # Baseline bacterial load set to one => all people are infected
            vaccinated=np.full(self.params['N'], fill_value=False, dtype=bool),
            time_since_vaccinated=np.zeros(self.params['N']),
            Infection_risk = np.random.gamma(size = self.params['N'], shape = self.params['infection_risk_shape'], scale = 1/self.params['infection_risk_shape'])
        )

    def getLambdaViaDirectCalculation(self, totalLoad, IndD, Infection_risk, vaccinated, time_since_vaccinated):
        prevLambda = self.bet * (self.params['v_1'] * totalLoad + self.params['v_2'] * (totalLoad ** (self.params['phi'] + 1)))
        a = 1/3
        b = 1/3
        c = 1/3
        epsm = 1 - self.params['epsilon']
        eps =  self.params['epsilon']
        returned = [eps * prevLambda[0] + prevLambda[0]*a * epsm + prevLambda[1]*epsm*b + prevLambda[2]*epsm*c,
                    prevLambda[0]*a*epsm + prevLambda[1]*b * epsm + eps * prevLambda[1] + prevLambda[2]*epsm*c,
                    prevLambda[0]*a*epsm + prevLambda[1]*epsm*b + prevLambda[2]*c * epsm + eps * prevLambda[2],
                ] * (0.5 + 0.5 * (1 - IndD)) * Infection_risk

        # add reduction in lambda according to who has been vaccinated
        prob_reduction = self.params["vacc_prob_block_transmission"]

        # add impact of waning using a linear slope. After waning period assumed vaccine has zero impact.
        prob_reduction = prob_reduction * (- time_since_vaccinated / self.params["vacc_waning_length"] + 1)
        prob_reduction = np.maximum(prob_reduction,0)

        returned[vaccinated] = (1 - prob_reduction[vaccinated]) * returned[vaccinated]


        return 1-np.exp(-returned)
    
    def test_WhenAllInfectedNoneDiseased(self):
        # First lambda_step calculation
        lambda_step = 1 - np.exp(-getlambdaStep(
            params=self.params, Age=self.Age, bact_load=self.vals['bact_load'],
            IndD=self.vals['IndD'], Infection_risk=self.vals['Infection_risk'], vaccinated=self.vals['vaccinated'],
            time_since_vaccinated=self.vals['time_since_vaccinated'], bet=self.bet, demog=self.demog
        ))
        # Assert array is as expected
        npt.assert_array_almost_equal(self.getLambdaViaDirectCalculation(totalLoad = np.array([1,1,1]),
                                                                         IndD = np.array([0,0,0]), 
                                                                         Infection_risk=self.vals['Infection_risk'],
                                                                         vaccinated = np.array([False,False,False]),
                                                                         time_since_vaccinated = np.array([0,0,0])),
                                    lambda_step)
    
    def test_WhenFirstPersonNotInfectedButIsDiseased(self):
        # Modify vals and recalculate lambda_step
        self.vals['bact_load'][0] = 0
        self.vals['IndD'][0] = 1

        lambda_step = 1 - np.exp(-getlambdaStep(
            params=self.params, Age=self.Age, bact_load=self.vals['bact_load'],
            IndD=self.vals['IndD'], Infection_risk=self.vals['Infection_risk'], vaccinated=self.vals['vaccinated'],
            time_since_vaccinated=self.vals['time_since_vaccinated'], bet=self.bet, demog=self.demog
        ))

        # Assert new values
        npt.assert_array_almost_equal(self.getLambdaViaDirectCalculation(totalLoad = np.array([0,1,1]),
                                                                         IndD = np.array([1,0,0]), 
                                                                         Infection_risk=self.vals['Infection_risk'],
                                                                         vaccinated = np.array([False,False,False]),
                                                                         time_since_vaccinated = np.array([0,0,0])),
                                    lambda_step)
        
    def test_WhenFirstPersonNotInfectedButIsDisease_SecondPersonWasVaccinatedLastWeek(self):
        # Change vaccinated and time_since_vaccinated and recalculate lambda_step
        self.vals['bact_load'][0] = 0
        self.vals['IndD'][0] = 1
        self.vals['time_since_vaccinated'] = np.array([0, 1, 0])
        self.vals['vaccinated'] = np.array([False, True, False])

        lambda_step = 1 - np.exp(-getlambdaStep(
            params=self.params, Age=self.Age, bact_load=self.vals['bact_load'],
            IndD=self.vals['IndD'], Infection_risk=self.vals['Infection_risk'], vaccinated=self.vals['vaccinated'],
            time_since_vaccinated=self.vals['time_since_vaccinated'], bet=self.bet, demog=self.demog
        ))

        # Assert new values
        npt.assert_array_almost_equal(self.getLambdaViaDirectCalculation(totalLoad = np.array([0,1,1]),
                                                                         IndD = np.array([1,0,0]), 
                                                                         Infection_risk=self.vals['Infection_risk'],
                                                                         vaccinated = np.array([False,True,False]),
                                                                         time_since_vaccinated = np.array([0,1,0])),
                                    lambda_step)
        
    def test_WhenFirstPersonNotInfectedButIsDisease_SecondPersonWasVaccinatedWaningLengthAgo(self):
        # Change vaccinated and time_since_vaccinated and recalculate lambda_step
        self.vals['bact_load'][0] = 0
        self.vals['IndD'][0] = 1
        # when the time_since_vaccinated = params["vacc_waning_length"], then it should be equivalent to
        # that person not being vaccinated. We will assume this when doing the direct calculation
        self.vals['time_since_vaccinated'] = np.array([0, self.params["vacc_waning_length"], 0])
        self.vals['vaccinated'] = np.array([False, True, False])
        
        lambda_step = 1 - np.exp(-getlambdaStep(
            params=self.params, Age=self.Age, bact_load=self.vals['bact_load'],
            IndD=self.vals['IndD'], Infection_risk=self.vals['Infection_risk'],
            vaccinated=self.vals['vaccinated'],
            time_since_vaccinated=self.vals['time_since_vaccinated'], bet=self.bet, demog=self.demog
        ))

        # Assert new values
        npt.assert_array_almost_equal(self.getLambdaViaDirectCalculation(totalLoad = np.array([0,1,1]),
                                                                         IndD = np.array([1,0,0]), 
                                                                         Infection_risk=self.vals['Infection_risk'],
                                                                         vaccinated = np.array([False,False,False]),
                                                                         time_since_vaccinated = np.array([0,self.params["vacc_waning_length"],0])),
                                    lambda_step)
