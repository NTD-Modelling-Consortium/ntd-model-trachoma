import numpy as np
from datetime import date
import pandas as pd
import pkg_resources


def readCoverageData(coverageFileName):
    # read coverage data file
    modelDataDir = pkg_resources.resource_filename( "trachoma", "data/coverage" )
    PlatCov = pd.read_csv(
         f"{modelDataDir}/{coverageFileName}"
    )
     # we want to find which is the first year specified in the coverage data, along with which
     # column of the data set this corresponds to
    fy_index = np.where(PlatCov.columns == "2020")[0][0]
    count = 0
    minAgeIndex = np.where(PlatCov.columns == "min age")[0][0]
    maxAgeIndex = np.where(PlatCov.columns == "max age")[0][0]
    for i in range(fy_index, len(PlatCov.columns)):
        dd = PlatCov.iloc[:, i]
        MDAS = np.where(dd>0)[0]
        if len(MDAS)>0:
            for k in range(len(MDAS)):
                j = MDAS[k]
                if count == 0:
                    MDAData = [[float(PlatCov.columns[i]), PlatCov.iloc[j, minAgeIndex], PlatCov.iloc[j, maxAgeIndex], PlatCov.iloc[j, i], j,  PlatCov.shape[0]]]
                    count += 1
                else:
                    MDAData.append([float(PlatCov.columns[i]), PlatCov.iloc[j, minAgeIndex], PlatCov.iloc[j, maxAgeIndex], PlatCov.iloc[j, i], j,  PlatCov.shape[0]])
                    count +=1
    return MDAData
                
                
                
def getMDADates(MDAData):
    for i in range(len(MDAData)):
        d = MDAData[i][0]
        y = int(d)
        m = round(12*(d - int(d))) + 1
        day = 1
        if i == 0:
            MDA_dates = [date(y, m, day)]
        else:
            MDA_dates.append(date(y, m, day))
    return MDA_dates



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



def Tx_matrix(params, sim_params, previous_rounds, MDAData = None):

    '''
    Create matrix to determine who gets treated at each MDA round,
    allowing for systematic non-compliance as specified by Dyson.
    If we include the MDAData input, then this will have specifications
    related to the MDAs and hence will be used to create the treatment 
    matrix. If this isn't included, we will use more basic information 
    related to treatment from the sim_params and params variables to
    create the treatment matrix.
    '''

    np.random.seed(0)

    if previous_rounds == 0:
        if MDAData is None:
            nMDA = sim_params['N_MDA']
            MDA_Cov = params['MDA_Cov']
            
        else:
            nMDA = len(MDAData)
            MDA_Cov = MDAData[0][3]
            
        # Assign first treatment
        ind_treat = np.zeros((params['N'], nMDA))
        ind_treat[:, 0] = np.random.uniform(size=params['N']) < MDA_Cov
        
        for k in range(1, nMDA):
            if MDAData is not None:
                MDA_Cov = MDAData[k][3]
            # Subsequent treatment probs function of previous treatments
            ind_treat[:, k] = np.random.binomial(n=1, size=params['N'], p=(MDA_Cov * (1 - params['rho']) +
            (params['rho'] * np.sum(ind_treat[:, :k], axis=1))) / (1 + (k + 1 - 2) * params['rho']))

    else:
        if MDAData is None:
            nMDA = sim_params['N_MDA']
            MDA_Cov = params['MDA_Cov']
        else:
            nMDA = len(MDAData)
            MDA_Cov = MDAData[0][3]
 
        ind_treat = np.zeros((params['N'], previous_rounds + nMDA))
        # Assign first treatment
        ind_treat[:, 0] = np.random.uniform(size=params['N']) < MDA_Cov

        for k in range(1, previous_rounds + nMDA):
            if MDAData is not None:
                MDA_Cov = MDAData[k][3]
            # Subsequent treatment probs function of previous treatments
            ind_treat[:, k] = np.random.binomial(n=1, size=params['N'], p=(MDA_Cov * (1 - params['rho']) +
            (params['rho'] * np.sum(ind_treat[:, :k], axis=1))) / (1 + (k + 1 - 2) * params['rho']))

        ind_treat = ind_treat[:, - nMDA:]

    return ind_treat



def get_MDA_times(MDA_dates, Start_date, burnin):
    MDA_times = []
    for i in range(0, len(MDA_dates)):
        MDA_times.append(burnin + int((MDA_dates[i] - Start_date).days/7))
    return np.array(MDA_times)



def getMDAAgeRanges(coverageFileName):

    modelDataDir = pkg_resources.resource_filename( "trachoma", "data/coverage" )
    PlatCov = pd.read_csv(
         f"{modelDataDir}/{coverageFileName}"
    )

    MDAAgeRanges = np.zeros([PlatCov.shape[0],2], dtype = object)
    minAgeIndex = np.where(PlatCov.columns == "min age")[0][0]
    maxAgeIndex = np.where(PlatCov.columns == "max age")[0][0]
    for i in range(PlatCov.shape[0]):
        MDAAgeRanges[i, 0] = PlatCov.iloc[i, minAgeIndex]
        MDAAgeRanges[i, 1] = PlatCov.iloc[i, maxAgeIndex]
        
    return MDAAgeRanges



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



def numMDAsBeforeNextSurvey(surveyPrev):
    '''
    Function to return the number of surveys before the next survey
    '''
    
    if surveyPrev >= 0.3:
        return 5 
    if surveyPrev >= 0.1:
        return 3
    if surveyPrev >= 0.05: 
        return 1
    return 100