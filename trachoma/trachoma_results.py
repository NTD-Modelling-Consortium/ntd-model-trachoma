#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import copy
from numpy import ndarray
from dataclasses import dataclass
from typing import Optional

@dataclass
class Result:
    time: float
    IndI: ndarray
    IndD: ndarray
    Age:ndarray
    NoInf: ndarray
    nMDA:Optional[ndarray] = None
    nMDADoses: Optional[ndarray] = None
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



def getResultsIHME(results, demog, params, outputYear):
    '''
    Function to collate results for IHME
    '''
    max_age = demog['max_age'] // 52 # max_age in weeks

    df = pd.DataFrame(0, range(len(outputYear)*4*60 ), columns= range(len(results)+4))
    df = df.rename(columns={0: "Time", 1: "age_start", 2: "age_end", 3: "measure"}) 
    
    for i in range(len(results)):
        ind = 0
        d = copy.deepcopy(results[i][1])
        for j in range(len(d)):
            year = outputYear[j]
            large_infection_count = (d[j].NoInf > params['n_inf_sev'])
            infection_count = (d[j].IndI > 0)
            Diseased = np.where(d[j].IndD == 1)
            NonDiseased = np.where(d[j].IndD == 0)
            pos = np.zeros(len(d[j].Age), dtype = object)
            if(len(Diseased) > 0):
                TruePositive = np.random.binomial(n=1, size=len(Diseased[0]), p = params['TestSensitivity'])
                pos[Diseased] = TruePositive
            if(len(NonDiseased) > 0):
                FalsePositive = np.random.binomial(n=1, size=len(NonDiseased[0]), p = 1- params['TestSpecificity'])
                pos[NonDiseased] = FalsePositive
            
 
            Age = d[j].Age
            # Cast weights to integer to be able to count
            manyInfs, _ = np.histogram(Age, bins=max_age, weights=large_infection_count.astype(int))
            Infs, _ = np.histogram(Age, bins=max_age, weights=infection_count.astype(int))
            nums, _ = np.histogram(Age, bins=max_age)
            observedDis, _ = np.histogram(Age, bins=max_age, weights=pos.astype(int))
            k = np.where(nums == 0)
            nums[k] = 1
            if i == 0:
                df.iloc[range(ind, ind+max_age), 0] = np.repeat(year,max_age)
                df.iloc[range(ind, ind+max_age), 1] = range(0, max_age)
                df.iloc[range(ind, ind+max_age), 2] = range(1, max_age + 1)
                df.iloc[range(ind, ind+max_age), 3] = np.repeat("TruePrevalence", max_age)
                df.iloc[range(ind, ind+max_age), i+4] = Infs/nums
                ind += max_age
                df.iloc[range(ind, ind+max_age), 0] = np.repeat(year,max_age)
                df.iloc[range(ind, ind+max_age), 1] = range(0, max_age)
                df.iloc[range(ind, ind+max_age), 2] = range(1, max_age + 1)
                df.iloc[range(ind, ind+max_age), 3] = np.repeat("ObservedTF", max_age)
                df.iloc[range(ind, ind+max_age), i+4] = observedDis/nums
                ind += max_age
                df.iloc[range(ind, ind+max_age), 0] = np.repeat(year,max_age)
                df.iloc[range(ind, ind+max_age), 1] = range(0, max_age)
                df.iloc[range(ind, ind+max_age), 2] = range(1, max_age + 1)
                df.iloc[range(ind, ind+max_age), 3] = np.repeat("heavyInfections", max_age)
                df.iloc[range(ind, ind+max_age), i+4] = manyInfs/nums
                ind += max_age
                nums[k] = 0
                df.iloc[range(ind, ind+max_age), 0] = np.repeat(year,max_age)
                df.iloc[range(ind, ind+max_age), 1] = range(0, max_age)
                df.iloc[range(ind, ind+max_age), 2] = range(1, max_age + 1)
                df.iloc[range(ind, ind+max_age), 3] = np.repeat("number", max_age)
                df.iloc[range(ind, ind+max_age), i+4] = nums
                ind += max_age
            else:
                df.iloc[range(ind, ind+max_age), i+4] = Infs/nums
                ind += max_age
                df.iloc[range(ind, ind+max_age), i+4] = observedDis/nums
                ind += max_age
                df.iloc[range(ind, ind+max_age), i+4] = manyInfs/nums
                ind += max_age
                df.iloc[range(ind, ind+max_age), i+4] = nums
                ind += max_age
    for i in range(len(results)):
        df = df.rename(columns={i+4: "draw_"+ str(i)}) 
    return df


def getResultsIPM(results, demog, params, outputYear, MDAAgeRanges):
    '''
    Function to collate results for IPM
    '''
   
    df = pd.DataFrame(0, range(len(outputYear)*3 + len(outputYear) * 3 * len(MDAAgeRanges)), columns= range(len(results)+4))
    df = df.rename(columns={0: "Time", 1: "age_start", 2: "age_end", 3: "measure"}) 
    
    for i in range(len(results)):
        ind = 0
        d = copy.deepcopy(results[i][1])
        for j in range(len(d)):
            year = outputYear[j]
            
            if i == 0:
                df.iloc[ind, 0] = year
                df.iloc[ind, 3] = "nSurvey"
                df.iloc[ind, 1] = "None"
                df.iloc[ind, 2] = "None"
                df.iloc[ind, i+4] = d[j].nSurvey
                ind += 1
                df.iloc[ind, 0] = year
                df.iloc[ind, 3] = "surveyPass"
                df.iloc[ind, 1] = "None"
                df.iloc[ind, 2] = "None"
                df.iloc[ind, i+4] = d[j].surveyPass
                ind += 1
                df.iloc[ind, 0] = year
                df.iloc[ind, 3] = "trueElimination"
                df.iloc[ind, 1] = "None"
                df.iloc[ind, 2] = "None"
                df.iloc[ind, i+4] = d[j].elimination
                ind += 1
                for k in range(len(MDAAgeRanges)):
                    df.iloc[ind, 0] = year
                    df.iloc[ind, 3] = "nDoses"
                    df.iloc[ind, 1] = MDAAgeRanges[k][0]
                    df.iloc[ind, 2] = MDAAgeRanges[k][1]
                    df.iloc[ind, i+4] = d[j].nMDADoses[k]
                    ind += 1
                    df.iloc[ind, 0] = year
                    df.iloc[ind, 3] = "MDAcoverage"
                    df.iloc[ind, 1] = MDAAgeRanges[k][0]
                    df.iloc[ind, 2] = MDAAgeRanges[k][1]
                    df.iloc[ind, i+4] = d[j].propMDA[k]
                    ind += 1
                    df.iloc[ind, 0] = year
                    df.iloc[ind, 3] = "numMDAs"
                    df.iloc[ind, 1] = MDAAgeRanges[k][0]
                    df.iloc[ind, 2] = MDAAgeRanges[k][1]
                    df.iloc[ind, i+4] = d[j].nMDA[k]
                    ind += 1
                    
                
            else:
                df.iloc[ind, i+4] = d[j].nSurvey
                ind += 1
                df.iloc[ind, i+4] = d[j].surveyPass
                ind += 1
                df.iloc[ind, i+4] = d[j].elimination
                ind += 1
                for k in range(len(MDAAgeRanges)):               
                    df.iloc[ind, i+4] = d[j].nMDADoses[k]
                    ind += 1   
                    df.iloc[ind, i+4] = d[j].propMDA[k]
                    ind += 1
                    df.iloc[ind, i+4] = d[j].nMDA[k]
                    ind += 1
                
    for i in range(len(results)):
        df = df.rename(columns={i+4: "draw_"+ str(i)}) 
    return df
