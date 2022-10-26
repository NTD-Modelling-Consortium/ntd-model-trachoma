import numpy as np
from datetime import date
import pandas as pd
import copy
import pkg_resources
from numpy import ndarray
from numpy.typing import NDArray
from dataclasses import dataclass
from typing import Callable, List, Optional

from trachoma.trachoma_results import *
from trachoma.trachoma_utils import *

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
    return vals



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



def doMDA(params, Age, MDA_round, Tx_mat, ageStart = None, ageEnd = None):

    '''
    Decide who is cured during MDA based on treatment matrix
    and probability of clearance given treated.
    '''
    if ageStart is None:
        babies = np.where(Age <= 26)[0]
        treated_babies = babies[Tx_mat[babies, MDA_round] == 1]
        cured_babies = treated_babies[np.random.uniform(size=len(treated_babies)) < (params['MDA_Eff'] * 0.5)]
    
        older = np.where(Age > 26)[0]
        treated_older = older[Tx_mat[older, MDA_round] == 1]
        cured_older = treated_older[np.random.uniform(size=len(treated_older)) < params['MDA_Eff']]
    else:
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
            treated_babies = []
    return np.append(cured_babies, cured_older), np.append(treated_babies, treated_older)



def MDA_timestep(vals, params, MDA_round, Tx_mat, ageStart = None, ageEnd = None):

    '''
    This is time step in which MDA occurs
    '''

    # Id who is treated and cured
    treated_cured = doMDA(params=params, Age=vals['Age'], MDA_round=MDA_round, Tx_mat=Tx_mat, ageStart = ageStart, ageEnd = ageEnd)

    # Set treated/cured indivs infection status and bacterial load to 0
    vals['IndI'][treated_cured[0].astype(int)] = 0       # clear infection they become I=0
    vals['bact_load'][treated_cured[0].astype(int)] = 0  # stop being infectious

    return vals, len(treated_cured[1])


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



def sim_Ind_MDA(params, Tx_mat, vals, timesim, demog, bet, MDA_times, MDAData = None, outputTimes = None, 
                               seed = 0, state = None, doSurvey = 0, outputRes = 0):

    '''
    Function to run a single simulation with MDA at time points determined by function MDA_times.
    Output is true prevalence of infection/disease in children aged 1-9.
    The Tx_mat variable contains information for each individual about if they will be included in treatment.
    In this function we have the option to run surveys of the population to decide if MDA should cease.
    This is governed by the doSurvey argument defaults to 0, and should be switched to 1 if we want to
    perform surveys.
    We can also include MDAData, which gives the specification of MDAs to be performed. If this is included,
    this will be used for age ranges of treatment. If not, then all individuals will be deemed as targets for 
    treatment.
    '''
    
    if outputTimes is not None:
        outputTimes2 = copy.deepcopy(outputTimes)    
        nextOutputTime = min(outputTimes2)
        w = np.where(outputTimes2 == nextOutputTime)
        outputTimes2[w] = timesim + 10
    else:
        nextOutputTime = timesim + 10
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
    yearly_threshold_infs = np.zeros((timesim + 1, int(demog['max_age']/52)))
   
   # get initial prevalence in 1-9 year olds. will decide how many MDAs (if any) to do before another survey
   
    if doSurvey == 1:
        surveyPrev = returnSurveyPrev(vals, params['TestSensitivity'], params['TestSpecificity'])
       # if the prevalence is <= 5%, then we have passed the survey and won't do any MDA
        #surveyPass = 0
        surveyPass = 1 if surveyPrev <= 0.05 else 0
       # if the prevalence is > 5%, then we will do another survey after given number of MDAs
       # call this value nextSurvey    
        nextSurvey = numMDAsBeforeNextSurvey(surveyPrev)
        # initialize time for next survey 
        surveyTime = min(MDA_times) + (nextSurvey * 52) + 26
        if surveyPrev <= 0.05:
            impactSurveyTime = min(MDA_times) + 104
        else:
            impactSurveyTime = timesim + 10
            
    #if MDAData is None:
        
   # initialize count of MDAs
    if MDAData is None:
        numMDA = np.zeros(1, dtype=object)
    else:
        numMDA = np.zeros(MDAData[0][-1], dtype=object)

    # initialize time for impact survey dependent on surveyed prevalence
    surveyPass = 0
    
    
    
    results = []
    nSurvey = 0
    prevNSurvey = 0
    if MDAData is None:
        nDoses = np.zeros(1, dtype=object)
        coverage = np.zeros(1, dtype=object)
        prevNMDA = np.zeros(1, dtype=object)
    else:
        nDoses = np.zeros(MDAData[0][-1], dtype=object)
        coverage = np.zeros(MDAData[0][-1], dtype=object)
        prevNMDA = np.zeros(MDAData[0][-1], dtype=object)
    
    for i in range(1, 1 + timesim):
        if i in MDA_times:
            if np.logical_or(surveyPass == 0, doSurvey == 0):
                MDA_round = np.where(MDA_times == i)[0]
                if(len(MDA_round) == 1):
                    MDA_round = MDA_round[0]
                    if MDAData is None:
                        ageStart = 0
                        ageEnd = 100
                    else:
                        ageStart = MDAData[MDA_round][1]
                        ageEnd = MDAData[MDA_round][2]
                    out = MDA_timestep(vals=vals, params=params, MDA_round=MDA_round, Tx_mat=Tx_mat, ageStart=ageStart, ageEnd=ageEnd)
                    vals = out[0]
                    if MDAData is None:
                        nDoses[0] += out[1]
                        # increment number of MDAs
                        numMDA[0] += 1
                        coverage[0] += out[1]/len(np.where(np.logical_and(vals['Age'] > ageStart * 52, vals['Age'] <= ageEnd *52))[0])
                        
                    else:
                        nDoses[MDAData[MDA_round][-2]] += out[1]
                        # increment number of MDAs
                        numMDA[MDAData[MDA_round][-2]] += 1
                        coverage[MDAData[MDA_round][-2]] += out[1]/len(np.where(np.logical_and(vals['Age'] > ageStart * 52, vals['Age'] <= ageEnd *52))[0])
                else:
                    for l in range(len(MDA_round)):
                        MDA_round2 = copy.deepcopy(MDA_round[l])
                        if MDAData is None:
                            ageStart = 0
                            ageEnd = 100
                        else:
                            ageStart = MDAData[MDA_round2][1]
                            ageEnd = MDAData[MDA_round2][2]
                        out = MDA_timestep(vals=vals, params=params, MDA_round=MDA_round2, Tx_mat=Tx_mat, ageStart=ageStart, ageEnd=ageEnd)
                        vals = out[0]
                        if MDAData is None:
                            nDoses[0] += out[1]
                            # increment number of MDAs
                            numMDA[0] += 1
                            coverage[0] += out[1]/len(np.where(np.logical_and(vals['Age'] > ageStart * 52, vals['Age'] <= ageEnd *52))[0])
                            
                        else:
                            nDoses[MDAData[MDA_round][-2]] += out[1]
                            # increment number of MDAs
                            numMDA[MDAData[MDA_round][-2]] += 1
                            coverage[MDAData[MDA_round][-2]] += out[1]/len(np.where(np.logical_and(vals['Age'] > ageStart * 52, vals['Age'] <= ageEnd *52))[0])
                # if the number of MDAs is the same as the number for the next survey then set survey time
          #  if sum(numMDA) == nextSurvey:
          #      surveyTime = i + 26
        #else:  removed and deleted one indent in the line below to correct mistake.
        #if np.logical_and(i == surveyTime, surveyPass==0):   
        if doSurvey == 1:
            if np.logical_or(i == surveyTime, i == impactSurveyTime) :     
                surveyPrev = returnSurveyPrev(vals, params['TestSensitivity'], params['TestSpecificity'])
                   
                # if the prevalence is <= 5%, then we have passed the survey and won't do any more MDA
                surveyPass = 1 if surveyPrev <= 0.05 else 0
                if surveyPass == 1:
                    impactSurveyTime = i + 104  
                # if the prevalence is > 5%, then we will do another survey after given number of MDAs
                # call this value nextSurvey    
                nextSurvey = numMDAsBeforeNextSurvey(surveyPrev)
                # add the number of MDAs already done to the number of MDAs to be done before the next survey
                surveyTime = i + (nextSurvey * 52) + 26
                
                nSurvey += 1
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
        # check if time to save variables to make Endgame outputs
        if np.logical_and(outputTimes is not None, i == nextOutputTime):
            
            # has the disease truly eliminated in the population
            true_elimination = 1 if (sum(vals['IndI']) + sum(vals['IndD'])) == 0 else 0
            # append the results to results variable
            results.append(outputResult(copy.deepcopy(vals), i, nDoses, coverage, numMDA-prevNMDA, nSurvey - prevNSurvey, surveyPass, true_elimination))
            # when will next Endgame output time be
            nextOutputTime = min(outputTimes2)
            # change next output time location in the all output variable to be after the end of the simulation
            # then the next output will be done at the correct time
            w = np.where(outputTimes2 == nextOutputTime)
            outputTimes2[w] = timesim + 10
            # save current num surveys, num MDAS as previous num surveys/MDAs, so next output we can tell how many were performed
            # since last output
            prevNSurvey = copy.deepcopy(nSurvey) 
            prevNMDA = copy.deepcopy(numMDA)
            # set coverage and nDoses to 0, so that if these are non-zero, we know that they occured since last output
            if MDAData is None:
                nDoses = np.zeros(1, dtype=object)
                coverage = np.zeros(1, dtype=object)
            else:    
                nDoses = np.zeros(MDAData[0][-1], dtype=object)
                coverage = np.zeros(MDAData[0][-1], dtype=object)
           
            
    vals['Yearly_threshold_infs'] = yearly_threshold_infs
    vals['True_Prev_Disease_children_1_9'] = prevalence # save the prevalence in children aged 1-9
    vals['True_Infections_Disease_children_1_9'] = infections # save the infections in children aged 1-9
    vals['State'] = np.random.get_state() # save the state of the simulations
    
    if outputRes == 0:
        return vals
    else:
        return vals, results







def run_single_simulation(pickleData, 
                          params, 
                          sim_params,
                          timesim,
                          demog, 
                          beta, 
                          MDA_times, 
                          MDAData, 
                          outputTimes, 
                          index,
                          doSurvey = 1,
                          outputRes = 1):
    
    '''
    Function to run a single instance of the simulation. The starting point for these simulations
    is
    '''
    vals = pickleData
    params['N'] = len(vals['IndI'])
    Tx_mat = Tx_matrix(params, sim_params, 0, MDAData)
    results = sim_Ind_MDA(params=params, Tx_mat = Tx_mat, 
                                        vals = vals, timesim = timesim,
                                        demog=demog, bet=beta, MDA_times = MDA_times, 
                                        MDAData=MDAData, outputTimes= outputTimes, 
                                        seed = index, doSurvey = doSurvey, 
                                        outputRes = outputRes)
    return results

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################




    
