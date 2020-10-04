import functools
import time
import numpy as np
import pandas as pd

from trachoma.trachoma_simulations import Trachoma_Simulation

def timer(func):
    """Print the runtime of the decorated function"""
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        print(f"-> Running {func.__name__!r}")
        start_time = time.perf_counter()    # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time    # 3
        print(f"=> Finished {func.__name__!r} in {run_time:.4f} secs\n\n")
        return value
    return wrapper_timer

@timer
def scenario_1():
    print("The simulations start in 2014 and end in 2020, while the MDA starts in 2015 and ends in 2019.")
    BetFilePath = 'files/InputBet_scen1.csv'  
    MDAFilePath = 'files/InputMDA_scen1.csv' 
    PrevFilePath = 'files/OutputPrev_scen1.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=False,
                        OutSimFilePath=None,
                        InSimFilePath=None)

    print("results in files/OutputPrev_scen1.csv")

@timer
def scenario_2():
    print("The simulations start in 2016 and end in 2020, while the MDA starts in 2017 and ends in 2019.")
    BetFilePath = 'files/InputBet_scen2.csv'  
    MDAFilePath = 'files/InputMDA_scen2.csv' 
    PrevFilePath = 'files/OutputPrev_scen2.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=False,
                        OutSimFilePath=None,
                        InSimFilePath=None)
        
    print("results in files/OutputPrev_scen2.csv")

@timer
def scenario_3():
    print("The simulations start in 2019 and end in 2020, while the there is no MDA.")
    BetFilePath = 'files/InputBet_scen3.csv'
    MDAFilePath = 'files/InputMDA_scen3.csv'
    PrevFilePath = 'files/OutputPrev_scen3.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=False,
                        OutSimFilePath=None,
                        InSimFilePath=None)
    print("results in files/OutputPrev_scen3.csv")

@timer
def scenario_4a():
    print( "The first set of simulations starts in 2018 and ends in 2020, while the MDA starts in 2019 and ends in 2020." )
    BetFilePath = 'files/InputBet_scen1.csv'
    MDAFilePath = 'files/InputMDA_scen4a.csv'
    PrevFilePath = 'files/OutputPrev_scen4a.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=True,
                        OutSimFilePath='files/OutputVals_scen4a.p',
                        InSimFilePath=None)

@timer
def scenario_4b():
    print( "The second set of simulations starts in 2021 and ends in 2023, while the MDA starts in 2021 and ends in 2022." )
    BetFilePath = 'files/InputBet_scen1.csv'
    MDAFilePath = 'files/InputMDA_scen4b.csv'
    PrevFilePath = 'files/OutputPrev_scen4b.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=False,
                        OutSimFilePath=None,
                        InSimFilePath='files/OutputVals_scen4a.p')

@timer
def scenario_4c():
    print( "For comparison purposes, the following code generates a third set of simulations starting in 2018 and ending in 2023, with MDA starting in 2019 and ending in 2022." )
    BetFilePath = 'files/InputBet_scen1.csv'  
    MDAFilePath = 'files/InputMDA_scen4c.csv' 
    PrevFilePath = 'files/OutputPrev_scen4c.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=False,
                        OutSimFilePath=None,
                        InSimFilePath=None)

@timer
def scenario_5a():
    print( "The first set of simulations starts in 2018 and ends in 2020, while there is no MDA." )
    BetFilePath = 'files/InputBet_scen1.csv'
    MDAFilePath = 'files/InputMDA_scen5a.csv'
    PrevFilePath = 'files/OutputPrev_scen5a.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=True,
                        OutSimFilePath='files/OutputVals_scen5a.p',
                        InSimFilePath=None)

@timer
def scenario_5b():
    print( "The second set of simulations starts in 2021 and ends in 2023, while there is still no MDA." )
    BetFilePath = 'files/InputBet_scen1.csv'  
    MDAFilePath = 'files/InputMDA_scen5b.csv' 
    PrevFilePath = 'files/OutputPrev_scen5b.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=False,
                        OutSimFilePath=None,
                        InSimFilePath='files/OutputVals_scen5a.p')

@timer
def scenario_5c():
    print( "For comparison purposes, the following code generates a third set of simulations starting in 2018 and ending in 2023, with no MDA." )
    BetFilePath = 'files/InputBet_scen1.csv'
    MDAFilePath = 'files/InputMDA_scen5c.csv'
    PrevFilePath = 'files/OutputPrev_scen5c.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=False,
                        OutSimFilePath=None,
                        InSimFilePath=None)

@timer
def scenario_6a():
    print( "The first set of simulations starts in 2018 and ends in 2020, while there is no MDA." )
    BetFilePath = 'files/InputBet_scen1.csv'  
    MDAFilePath = 'files/InputMDA_scen6a.csv' 
    PrevFilePath = 'files/OutputPrev_scen6a.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=True,
                        OutSimFilePath='files/OutputVals_scen6a.p',
                        InSimFilePath=None)

@timer
def scenario_6b():
    print( "The second set of simulations starts in 2021 and ends in 2023, with MDA also from 2021 to 2023." )
    BetFilePath = 'files/InputBet_scen1.csv'
    MDAFilePath = 'files/InputMDA_scen6b.csv'
    PrevFilePath = 'files/OutputPrev_scen6b.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=False,
                        OutSimFilePath=None,
                        InSimFilePath='files/OutputVals_scen6a.p')

@timer
def scenario_6c():
    print( "For comparison purposes, the following code generates a third set of simulations starting in 2018 and ending in 2023, with MDA from 2021 to 2023." )
    BetFilePath = 'files/InputBet_scen1.csv'
    MDAFilePath = 'files/InputMDA_scen6c.csv'
    PrevFilePath = 'files/OutputPrev_scen6c.csv'

    Trachoma_Simulation(BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        SaveOutput=False,
                        OutSimFilePath=None,
                        InSimFilePath=None)

scenario_1()
scenario_2()
scenario_3()
scenario_4a()
scenario_4b()
scenario_4c()
scenario_5a()
scenario_5b()
scenario_5c()
scenario_6a()
scenario_6b()
scenario_6c()
