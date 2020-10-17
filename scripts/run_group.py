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
                'percentile_25': op_data.iloc[:, 2:].quantile(0.25),
                'percentile_75': op_data.iloc[:, 2:].quantile(0.75)
            }).to_json( summary_json_path )

print( f"===== FINISHED RUNNING GROUP {group} =====" )
