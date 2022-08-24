import numpy as np
from datetime import date
import matplotlib.pyplot as plt
import pandas as pd
import copy
import pkg_resources
from numpy import ndarray
from numpy.typing import NDArray
from dataclasses import dataclass
from typing import Callable, List, Optional


@dataclass
class Result:
    time: float
    IndI: ndarray
    IndD: ndarray
    Age:ndarray
    NoInf: ndarray
    nMDA:Optional[int] = None
    nMDADoses: Optional[int] = None
    nSurvey: Optional[int] = None
    surveyPass: Optional[int] = None
    elimination: Optional[int] = None
    propMDA: Optional[ndarray] = None
    
def outputResult(vals, i, nDoses, coverage, nMDA, nSurvey, surveyPass, true_elimination):
    return (Result(time = i,
                          IndI = vals['IndI'], 
                          IndD = vals['IndD'], 
                          Age = vals['Age'], 
                          NoInf = vals['No_Inf'],
                          nMDADoses = nDoses, 
                          nSurvey = nSurvey,
                          surveyPass = surveyPass,
                          elimination = true_elimination,
                          propMDA = coverage,
                          nMDA = nMDA))

def stepF_fixed(vals, params, demog, bet):

    '''
    Step function i.e. transitions in each time non-MDA timestep.
    '''
    # Step 1: Identify individuals available for infection.
    # Susceptible individuals available for infection.
    Ss = np.where(vals['IndI'] == 0)[0]

    # Step 2: Calculate infection pressure from previous time step and choose infected individuals
    # Susceptible individuals acquiring new infections. This gives a lambda
    # for each individual dependent on age and disease status.
    lambda_step = 1 - np.exp(- getlambdaStep(params=params, Age=vals['Age'], bact_load=vals['bact_load'],
    IndD=vals['IndD'], bet=bet, demog=demog))
    # New infections
    newInf = Ss[np.random.uniform(size=len(Ss)) < lambda_step[Ss]]

    # Step 3: Identify transitions
    newDis = np.where(vals['T_latent'] == 1)[0]  # Designated latent period for that individual is about to expire
    newClearInf = np.where(vals['T_ID'] == 1)[0]  # Designated infectious period for that individual is about to expire
    newClearDis = np.where(vals['T_D'] == 1)[0]  # Designated diseased period for that individual is about to expire
    newInfectious = np.where(np.logical_and(vals['IndI']==1,vals['T_latent']==1))[0] # Only individuals who have avoided MDA become infectious at end of latent

    # Step 4: reduce counters
    # Those in latent period should count down with each timestep.
    vals['T_latent'][vals['T_latent'] > 0] -= 1
    # Those infected should count down with each timestep
    vals['T_ID'][vals['T_ID'] > 0] -= 1
    # Those diseased should count down with each timestep
    vals['T_D'][vals['T_D'] > 0] -= 1

    # Step 5: implement transitions
    # Transition: become diseased (and infected)
    vals['IndD'][newDis] = 1  # if they've become diseased they become D=1
    vals['T_ID'][newDis] = ID_period_function(Ind_ID_period_base=vals['Ind_ID_period_base'][newDis],
    No_Inf=vals['No_Inf'][newDis], params=params)
    #vals['T_D'][newDis] = 0  # SS Added to prevent transition of doom.
    # Transition: Clear infection
    vals['IndI'][newClearInf] = 0  # clear infection they become I=0
    # When individual clears infection, their diseased only is set
    vals['T_D'][newClearInf] = D_period_function(Ind_D_period_base=vals['Ind_D_period_base'][newClearInf],
    No_Inf=vals['No_Inf'][newClearInf], params=params, Age = vals['Age'][newClearInf])
    # Stop being infectious too
    vals['bact_load'][newClearInf] = 0
    # Transition: Clear disease
    vals['IndD'][newClearDis] = 0  # clear disease they become D=0
    # Transition: Become infectious
    vals['bact_load'][newInfectious] = bacterialLoad(No_Inf=vals['No_Inf'][newInfectious])

    # Step 6: implement infections
    # Transition: become infected
    vals['IndI'][newInf] = 1  # if they've become infected, become I=1
    # When individual becomes infected, set their latent period;
    # this is how long they remain in category I (infected but not diseased)
    vals['T_latent'][newInf] = vals['Ind_latent'][newInf]
    # New infected can be D and have nonzero T_D/T_ID
    vals['T_D'][newInf] = 0
    vals['T_ID'][newInf] = 0

    # Tracking infection history
    vals['No_Inf'][newInf] += 1

    # Update age, all age by 1w at each timestep, and resetting all "reset indivs" age to zero
    # Reset_indivs - Identify individuals who die in this timestep, either reach max age or random death rate
    vals['Age'] += 1
    reset_indivs = Reset(Age=vals['Age'], demog=demog, params=params)

    # Resetting new parameters for all new individuals created
    vals['Age'][reset_indivs] = 0
    vals['IndI'][reset_indivs] = 0
    vals['IndD'][reset_indivs] = 0
    vals['No_Inf'][reset_indivs] = 0
    vals['T_latent'][reset_indivs] = 0
    vals['T_ID'][reset_indivs] = 0
    vals['T_D'][reset_indivs] = 0
    #me = 2
    #print(vals['Age'][me],vals['No_Inf'][me],vals['bact_load'][me],':',vals['IndI'][me],vals['IndD'][me],vals['T_latent'][me],vals['T_ID'][me],vals['T_D'][me])

    return vals



def get_MDA_times(MDA_dates, Start_date, burnin):
    MDA_times = []
    for i in range(0, len(MDA_dates)):
        MDA_times.append(burnin + int((MDA_dates[i] - Start_date).days/7))
    return np.array(MDA_times)


def getOutputTimes(outputTimes):
    for i in range(len(outputTimes)):
        d = outputTimes[i]
        y = int(d)
        m = 6
        day = 1
        if i == 0:
            modOutputTimes = [date(y, m, day)]
        else:
            modOutputTimes.append(date(y, m, day))
    return modOutputTimes

def getlambdaStep(params, Age, bact_load, IndD, bet, demog):

    y_children = np.where(np.logical_and(Age >= 0, Age < 9 * 52))[0]  # Young children
    o_children = np.where(np.logical_and(Age >= 9 * 52, Age < 15 * 52))[0]  # Older children
    adults = np.where(Age >= 15 * 52)[0]  # Adults

    totalLoad = np.array([np.sum(bact_load[y_children]) / len(y_children),
    np.sum(bact_load[o_children]) / len(o_children), np.sum(bact_load[adults]) / len(adults)])
    prevLambda = bet * (params['v_1'] * totalLoad + params['v_2'] * (totalLoad ** (params['phi'] + 1)))

    a = len(y_children)/params['N']
    b = len(o_children)/params['N']
    c = len(adults)/params['N']
    epsm = 1 - params['epsilon']
    A = [
        prevLambda[0]*a + prevLambda[1]*epsm*b + prevLambda[2]*epsm*c,
        prevLambda[0]*a*epsm + prevLambda[1]*b + prevLambda[2]*epsm*c,
        prevLambda[0]*a*epsm + prevLambda[1]*epsm*b + prevLambda[2]*c,
    ]
    returned = np.ones(params['N'])
    returned[y_children] = A[0]
    returned[o_children] = A[1]
    returned[adults] = A[2]
    return returned

def Reset(Age, demog, params):

    '''
    Function to identify individuals who either die due
    to background mortality, or who reach max age.
    '''

    return np.where(np.logical_or(np.random.uniform(size=params['N']) < 1 - np.exp(- demog['tau']), Age > demog['max_age']))[0]

def Tx_matrix(params, sim_params, previous_rounds):

    '''
    Create matrix to determine who gets treated at each MDA round,
    allowing for systematic non-compliance as specified by Dyson.
    '''

    np.random.seed(0)

    if previous_rounds == 0:

        ind_treat = np.zeros((params['N'], sim_params['N_MDA']))

        # Assign first treatment
        ind_treat[:, 0] = np.random.uniform(size=params['N']) < params['MDA_Cov']

        for k in range(1, sim_params['N_MDA']):

            # Subsequent treatment probs function of previous treatments
            ind_treat[:, k] = np.random.binomial(n=1, size=params['N'], p=(params['MDA_Cov'] * (1 - params['rho']) +
            (params['rho'] * np.sum(ind_treat[:, :k], axis=1))) / (1 + (k + 1 - 2) * params['rho']))

    else:

        ind_treat = np.zeros((params['N'], previous_rounds + sim_params['N_MDA']))

        # Assign first treatment
        ind_treat[:, 0] = np.random.uniform(size=params['N']) < params['MDA_Cov']

        for k in range(1, previous_rounds + sim_params['N_MDA']):

            # Subsequent treatment probs function of previous treatments
            ind_treat[:, k] = np.random.binomial(n=1, size=params['N'], p=(params['MDA_Cov'] * (1 - params['rho']) +
            (params['rho'] * np.sum(ind_treat[:, :k], axis=1))) / (1 + (k + 1 - 2) * params['rho']))

        ind_treat = ind_treat[:, - sim_params['N_MDA']:]

    return ind_treat

def doMDA(params, Age, MDA_round, Tx_mat):

    '''
    Decide who is cured during MDA based on treatment matrix
    and probability of clearance given treated.
    '''

    babies = np.where(Age <= 26)[0]
    treated_babies = babies[Tx_mat[babies, MDA_round] == 1]
    cured_babies = treated_babies[np.random.uniform(size=len(treated_babies)) < (params['MDA_Eff'] * 0.5)]

    older = np.where(Age > 26)[0]
    treated_older = older[Tx_mat[older, MDA_round] == 1]
    cured_older = treated_older[np.random.uniform(size=len(treated_older)) < params['MDA_Eff']]
    #print('MDA:',cured_babies,cured_older)
    return np.append(cured_babies, cured_older)


def doMDAAgeRange(params, Age, MDA_round, Tx_mat, ageStart, ageEnd):

    '''
    Decide who is cured during MDA based on treatment matrix
    and probability of clearance given treated.
    '''
    if ageStart*52 <= 26:
        babies = np.where(Age <= 26)[0]
        treated_babies = babies[Tx_mat[babies, MDA_round] == 1]
        cured_babies = treated_babies[np.random.uniform(size=len(treated_babies)) < (params['MDA_Eff'] * 0.5)]
    
        older = np.where(np.logical_and(Age > 26, Age <= ageEnd *52))[0]
        treated_older = older[Tx_mat[older, MDA_round] == 1]
        cured_older = treated_older[np.random.uniform(size=len(treated_older)) < params['MDA_Eff']]
    else:
        older = np.where(np.logical_and(Age > ageStart * 52, Age <= ageEnd *52))[0]
        treated_older = older[Tx_mat[older, MDA_round] == 1]
        cured_older = treated_older[np.random.uniform(size=len(treated_older)) < params['MDA_Eff']]
        cured_babies = []
    #print('MDA:',cured_babies,cured_older)
    return np.append(cured_babies, cured_older)

def MDA_timestep(vals, params, MDA_round, Tx_mat):

    '''
    This is time step in which MDA occurs
    '''

    # Id who is treated and cured
    treated_cured = doMDA(params=params, Age=vals['Age'], MDA_round=MDA_round, Tx_mat=Tx_mat)

    # Set treated/cured indivs infection status and bacterial load to 0
    vals['IndI'][treated_cured] = 0       # clear infection they become I=0
    vals['bact_load'][treated_cured] = 0  # stop being infectious

    return vals, len(treated_cured)

def MDA_timestep_Age_range(vals, params, MDA_round, Tx_mat, ageStart, ageEnd):

    '''
    This is time step in which MDA occurs
    '''

    # Id who is treated and cured
    treated_cured = doMDAAgeRange(params=params, Age=vals['Age'], MDA_round=MDA_round, Tx_mat=Tx_mat, ageStart = ageStart, ageEnd = ageEnd)

    # Set treated/cured indivs infection status and bacterial load to 0
    vals['IndI'][treated_cured] = 0       # clear infection they become I=0
    vals['bact_load'][treated_cured] = 0  # stop being infectious

    return vals, len(treated_cured)


def ID_period_function(Ind_ID_period_base, No_Inf, params):

    '''
    Function to give duration of active infection.
    '''

    return np.round((Ind_ID_period_base - params['min_ID']) * np.exp(-params['inf_red'] * (No_Inf - 1)) + params['min_ID'])

def D_period_function(Ind_D_period_base, No_Inf, params, Age):

    '''
    Function to give duration of disease only period.
    '''
    ag = 0.00179
    aq = 0.0368
    T_ID = np.round( params['min_D'] + ( Ind_D_period_base - params['min_D'] ) * np.exp( -aq * ( No_Inf-1 ) -ag * Age ) )

    return T_ID

def bacterialLoad(No_Inf):

    '''
    Function to scale bacterial load according to infection history
    '''

    b1 = 1
    ep2 = 0.114

    return b1 * np.exp((No_Inf - 1) * - ep2)

def Set_inits(params, demog, sim_params):

    '''
    Set initial values.
    '''

    np.random.seed(0)

    vals = dict(

        # Individual's infected status
        IndI=np.zeros(params['N']),

        # Individual's disease status
        IndD=np.zeros(params['N']),

        # Individual's total number of infections, should increase by 1 each time they become newly infected
        No_Inf=np.zeros(params['N']),

        # Duration of latent period (I), i.e. infected but not yet clinically diseased
        T_latent=np.zeros(params['N']),

        # Duration of current ID period, set when becomes infected and counts down with each time step
        T_ID=np.zeros(params['N']),

        # Duration individual spends diseased after clearing infection
        T_D=np.zeros(params['N']),

        # Individual's latent period (fixed for now)
        Ind_latent=params['av_I_duration'] * np.ones(params['N']),

        # Individual's baseline ID period (first infection)
        Ind_ID_period_base=np.random.poisson(lam=params['av_ID_duration'], size=params['N']),

        # Individual's baseline diseased period (first infection)
        Ind_D_period_base=np.random.poisson(lam=params['av_D_duration'], size=params['N']),

        # Baseline bacterial load set to zero
        bact_load=np.zeros(params['N']),

        # Age distribution
        Age=init_ages(params=params, demog=demog),

        # Number of MDA rounds
        N_MDA=sim_params['N_MDA'],

        # Prevalence
        True_Prev_Disease_children_1_9=[],

    )

    return vals

def Seed_infection(params, vals):

    '''
    Seed infection.
    '''

    # set 1% to infected, can change if want to seed more, may need to if want to simulate
    # low transmission settings to stop stochastic fade-out during burn-in
    vals['IndI'][:int(np.round(params['N'] * 0.01))] = 1

    Init_infected = np.where(vals['IndI'] == 1)[0]

    # set latent period for those infected at start of simulation
    vals['T_latent'][Init_infected] = vals['Ind_latent'][Init_infected]

    # set number of infections to 1 for those infected at start of simulation
    vals['No_Inf'][Init_infected] = 1

    return vals

def init_ages(params, demog):

    '''
    Initialise age distribution
    Note: ages are in weeks.
    '''

    np.random.seed(0)

    ages = np.arange(1, 1 + demog['max_age'])

    # ensure the population is in equilibrium
    propAges = np.empty(len(ages))
    propAges[:-1] = np.exp(-ages[:-1] / demog['mean_age']) - np.exp(-ages[1:] / demog['mean_age'])
    propAges[-1] = 1 - np.sum(propAges[:-1])

    return np.random.choice(a=ages, size=params['N'], replace=True, p=propAges)

def sim_Ind_MDA(params, Tx_mat, vals, timesim, demog, bet, MDA_times, seed, state=None):

    '''
    Function to run a single simulation with MDA at time points determined by function MDA_times.
    Output is true prevalence of infection/disease in children aged 1-9.
    '''

    # when we are starting new simulations
    # we use the provided random seed
    if state is None:

        np.random.seed(seed)

    # when we are resuming previous simulations
    # we use the provided random state
    else:

        np.random.set_state(state)

    prevalence = []
    infections = []
    max_age = demog['max_age'] // 52 # max_age in weeks
    yearly_threshold_infs = np.zeros((timesim+1, max_age))
    for i in range(1, 1 + timesim):

        if i in MDA_times:

            MDA_round = np.where(MDA_times == i)[0][0]

            out = MDA_timestep(vals=vals, params=params, MDA_round=MDA_round, Tx_mat=Tx_mat)
            vals = out[0]
            nDoses = out[1]
            coverage = nDoses/ len(vals['IndI'])

        #else:  removed and deleted one indent in the line below to correct mistake.

        vals = stepF_fixed(vals=vals, params=params, demog=demog, bet=bet)

        children_ages_1_9 = np.logical_and(vals['Age'] < 10 * 52, vals['Age'] >= 52)
        n_children_ages_1_9 = np.count_nonzero(children_ages_1_9)
        n_true_diseased_children_1_9 = np.count_nonzero(vals['IndD'][children_ages_1_9])
        n_true_infected_children_1_9 = np.count_nonzero(vals['IndI'][children_ages_1_9])
        prevalence.append(n_true_diseased_children_1_9 / n_children_ages_1_9)
        infections.append(n_true_infected_children_1_9 / n_children_ages_1_9)

        large_infection_count = (vals['No_Inf'] > params['n_inf_sev'])
        # Cast weights to integer to be able to count
        a, _ = np.histogram(vals['Age'], bins=max_age, weights=large_infection_count.astype(int))
        yearly_threshold_infs[i, :] = a / params['N']

    vals['Yearly_threshold_infs'] = yearly_threshold_infs
    vals['True_Prev_Disease_children_1_9'] = prevalence # save the prevalence in children aged 1-9
    vals['True_Infections_Disease_children_1_9'] = infections # save the infections in children aged 1-9
    vals['State'] = np.random.get_state() # save the state of the simulations

    return vals


def numMDAsBeforeNextSurvey(surveyPrev):
    '''
    Function to return the number of MDAs before the next survey. Based on WHO guidelines
    dependent on the prevalence in the population, we will do a certain number of MDAs
    before checking the prevalence again
    '''
    
    if surveyPrev >= 0.3:
        nextSurvey = 5
    elif surveyPrev >= 0.1 & surveyPrev < 0.3:
        nextSurvey = 3
    elif surveyPrev >= 0.05 & surveyPrev < 0.1:
        nextSurvey = 1
    else:
        nextSurvey = 0
        
    return nextSurvey

def sim_Ind_MDA_Include_Survey(params, Tx_mat, vals, timesim, demog, bet, MDA_times, seed, state=None):

    '''
    Function to run a single simulation with MDA at time points determined by function MDA_times.
    Output is true prevalence of infection/disease in children aged 1-9.
    '''

    # when we are starting new simulations
    # we use the provided random seed
    if state is None:

        np.random.seed(seed)

    # when we are resuming previous simulations
    # we use the provided random state
    else:

        np.random.set_state(state)

    prevalence = []
    infections = []
    yearly_threshold_infs = np.zeros(( timesim+1, int(demog['max_age']/52)))
    for i in range(1, 1 + timesim):

        if i in MDA_times:

            MDA_round = np.where(MDA_times == i)[0][0]

            out = MDA_timestep(vals=vals, params=params, MDA_round=MDA_round, Tx_mat=Tx_mat)
            vals = out[0]
            nDoses = out[1]
            coverage = nDoses/ len(vals['IndI'])

        #else:  removed and deleted one indent in the line below to correct mistake.

        vals = stepF_fixed(vals=vals, params=params, demog=demog, bet=bet)

        a = []
        children_ages_1_9 = np.logical_and(vals['Age'] < 10 * 52, vals['Age'] >= 52)
        n_children_ages_1_9 = children_ages_1_9.sum()
        n_true_diseased_children_1_9 = vals['IndD'][children_ages_1_9].sum()
        n_true_infected_children_1_9 = vals['IndI'][children_ages_1_9].sum()
        prevalence.append(n_true_diseased_children_1_9 / n_children_ages_1_9)
        infections.append(n_true_infected_children_1_9 / n_children_ages_1_9)
        for l in range(0, int(demog['max_age']/52)):
            index = np.logical_and(vals['Age'] >= (l*52), vals['Age'] < ((l+1)*52))
            p = sum(vals['No_Inf'][index] > params['n_inf_sev'])/len(index)
            a.append(p)

        yearly_threshold_infs[i, :] = a

    vals['Yearly_threshold_infs'] = yearly_threshold_infs
    vals['True_Prev_Disease_children_1_9'] = prevalence # save the prevalence in children aged 1-9
    vals['True_Infections_Disease_children_1_9'] = infections # save the infections in children aged 1-9
    vals['State'] = np.random.get_state() # save the state of the simulations

    return vals





def returnSurveyPrev(vals, TestSensitivity, TestSpecificity):
    '''
    Function to run a return the tested prevalence of 1-9 year olds.
    This includes sensitivity and specificity of the test.
    Will be used in surveying to decide if we should do MDA, and how many MDAs before next test
    '''
    # survey 1-9 year olds
    children_ages_1_9 = np.logical_and(vals['Age'] < 10 * 52, vals['Age'] >= 52)
    
    # calculate true number of diseased 1-9 year olds
    Diseased = vals['IndD'][children_ages_1_9].sum()
    
    # how many 1-9 year olds are not diseased
    NonDiseased = children_ages_1_9.sum() - Diseased
    
    # perform test with given sensitivity and specificity to get test positives
    positive = int(np.random.binomial(n=Diseased, size=1, p = TestSensitivity)) + int(np.random.binomial(n=NonDiseased, size=1, p = 1- TestSpecificity)) 
    
    # return prevalence,calculated as number who tested positive divided by number of 1-9 year olds
    return positive/children_ages_1_9.sum()

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
