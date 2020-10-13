'''
Panayiota has grouped the IUs into groups so we need to do less runs. To run each group we would need to:
- take the InputBet_groupx file which contains the 200 beta parameters and pass it as BetFilePath
- take the OutputVals_groupx file which contains the population and pass it as InSimFilePath
- the PrevFilePath is just where the model outputs the runs so this should be where we define the naming etc for the output
- SaveOutput can be false and OutSimFilePath can be none

This is the basic setup for a group of ius and the foundation for our various scenarios. For the scenarios we want to have results for:
- MDAs on or off every 6 months for 10 years (2020 - 2030)
- and we would need to run all of these scenarios for rho 0, 0.1, 0.2, 0.3 
(Rho is currently a fixed value in https://github.com/ArtRabbitStudio/ntd-model-trachoma/blob/master/trachoma/trachoma_simulations.py on line 161, we need to write something so we can pass this in when we run our scenarios)

These are the group scenarios we need, we then have the mapping from groups to actual IUs which is in IUcodeTrachomaGroupsEthiopia.csv in Panaiyotas files. The relevant columns are IUcode and Group.
That’s all the scenarios we need for the simulator output.

For the setup page where I want to display historical prevalence data we’d need to take the OutputPrev_groupx file take the prevalence for each year (eg. 02-2014) and find the mean of the 200 runs. I just need a row with prevalence values and have all the IUs added together into the country-level.csv file (example in public/diseases/trachoma/country-level.csv).

We should also think about that she’s just given us Ethiopia and can give us the other countries as well, so what we write should be able to take another set of files, or once we got it going and everything seems to be good, we ask her to produce a full set.
'''

import itertools
import pandas as pd
import json
import pathlib
import sys
import os
from trachoma.trachoma_simulations import Trachoma_Simulation

if len( sys.argv ) != 2 or sys.argv[ 1 ].isdigit() == False:
    print( f"usage: {sys.argv[0]} <group number>" )
    sys.exit()

group = sys.argv[ 1 ]

annual = [ s for s in range( 202001, 203101, 100 ) ]
biannual = sorted( annual + [ s for s in range( 202006, 203106, 100 ) ] )

mda_vectors = {
    '12': annual,
    '6': biannual
}

FilePathRoot = f"data/merged"
BetFilePath = f"{FilePathRoot}/InputBet_group{group}.csv"  
InSimFilePath = f"{FilePathRoot}/OutputVals_group{group}.p"  

if not pathlib.Path( BetFilePath ).is_file():
    print( f"Input beta file {BetFilePath} doesn't exist, exiting." )
    sys.exit()

if not pathlib.Path( InSimFilePath ).is_file():
    print( f"Input simulation file {InSimFilePath} doesn't exist, exiting." )
    sys.exit()

# specified MDA coverage values
for MDA_Cov in [ 0.6, 0.7, 0.8, 0.9 ]:

    # annual, biannual
    for mda_type, mda_vector in mda_vectors.items():

        mda_lists = [ mda_vector[i:j] for i, j in itertools.combinations( range( len( mda_vector ) + 1 ), 2 ) ]

        for mda_list in mda_lists:

            # mkdir -p group dirs
            output_dir = f"{FilePathRoot}/output/group-{group}/coverage-{MDA_Cov}/mdatype-{mda_type}"
            pathlib.Path( output_dir ).mkdir( parents = True, exist_ok = True )

            input_dir = f"{FilePathRoot}/input/group-{group}/coverage-{MDA_Cov}/mdatype-{mda_type}"
            pathlib.Path( input_dir ).mkdir( parents = True, exist_ok = True )

            mda_list_string = '-'.join( [ str( x ) for x in mda_list ] )
            file_name_root = f"{group}-{MDA_Cov}-{mda_type}-{mda_list_string}"

            # Input MDA file path/name
            input_mda_file_name = f"InputMDA-{file_name_root}.csv"
            MDAFilePath = f"{input_dir}/{input_mda_file_name}"

            # Input MDA data
            df = pd.DataFrame.from_records( [ {
                'start_sim_year': 2020,
                'end_sim_year': 2030,
                'first_mda': '',
                'last_mda': '',
                'mda_vector': json.dumps( mda_list )
            } ] )

            # write Input MDA to file
            df.to_csv( MDAFilePath, index=None )

            # output CSV file path
            csv_file_name = f"{file_name_root}.csv"
            PrevFilePath = f"{output_dir}/{csv_file_name}"

            print( f"=== Running:\n\tMDA_Cov {MDA_Cov}\n\tvector {mda_type}\n\tinput {MDAFilePath}\n\toutput {PrevFilePath}\n" )
            continue

            # run the simulation
            Trachoma_Simulation(
                BetFilePath=BetFilePath,
                MDAFilePath=MDAFilePath,
                PrevFilePath=PrevFilePath,
                SaveOutput=False,
                OutSimFilePath=None,
                InSimFilePath=InSimFilePath,
                MDA_Cov=MDA_Cov
            )

            # remove the input file
            os.remove( MDAFilePath )

            # read in the simulation output
            op_data = pd.read_csv( PrevFilePath )

            # make a json file path to summarise it into
            json_file_name = f"{file_name_root}-summary.json"
            summary_json_path = f"{output_dir}/{json_file_name}"

            # summarise it in there 
            pd.DataFrame( {
                'median': op_data.iloc[:, 2:].median(),
                '25_percentile': op_data.iloc[:, 2:].quantile(0.25),
                '75_percentile': op_data.iloc[:, 2:].quantile(0.75)
            }).to_json( summary_json_path )
