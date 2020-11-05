import itertools
import pandas as pd
import json
import pathlib
import sys
import os
import re
from google.cloud import storage
from trachoma.trachoma_simulations import Trachoma_Simulation

if len( sys.argv ) != 2 :
    print( f"usage: {sys.argv[0]} <group number>" )
    sys.exit()

groups = sys.argv[ 1 ].split( ',' )

annual = [ s for s in range( 202001, 203101, 100 ) ]
biannual = sorted( annual + [ s for s in range( 202006, 203106, 100 ) ] )

mda_vectors = {
    '12': annual,
    '6': biannual
}

FilePathRoot = f"data/merged"
CloudPathRoot = f"diseases/trachoma/data"

client = storage.Client()
bucket = client.get_bucket('ntd-disease-simulator-data')

for group in groups:

    ''' SET UP INPUT FILE PATHS '''

    BetFilePath = f"{FilePathRoot}/InputBet_group{group}.csv"
    InSimFilePath = f"{FilePathRoot}/OutputVals_group{group}.p"

    if not pathlib.Path( BetFilePath ).is_file():
        print( f"Input beta file {BetFilePath} doesn't exist, exiting." )
        sys.exit()

    if not pathlib.Path( InSimFilePath ).is_file():
        print( f"Input simulation file {InSimFilePath} doesn't exist, exiting." )
        sys.exit()

    GroupOutputFileDir = f"{FilePathRoot}/output/group-{group}"
    pathlib.Path( GroupOutputFileDir ).mkdir( parents = True, exist_ok = True )

    ''' SUMMARISE HISTORICAL DATA '''

    # read in historical data file
    HistPrevFilePath = f"{FilePathRoot}/OutputPrev_group{group}.csv"
    historical_prevalence = pd.read_csv( HistPrevFilePath )

    # make a json file path to summarise it into
    hist_summary_file_name = f"{group}-historical-prevalence-summary.json"
    HistSummaryFilePath = f"{GroupOutputFileDir}/{hist_summary_file_name}"

    # summarise it in there
    pd.DataFrame( {
        'median': historical_prevalence.iloc[:, 2:].median(),
        'lower': historical_prevalence.iloc[:, 2:].quantile(0.05),
        'upper': historical_prevalence.iloc[:, 2:].quantile(0.95)
    }).to_json( HistSummaryFilePath )

    ''' UPLOAD HISTORICAL DATA & SUMMARY '''

    # upload summary file
    HistSummaryCloudPath = f"{CloudPathRoot}/group-{group}/{hist_summary_file_name}"
    historical_summary_blob = bucket.blob( HistSummaryCloudPath )
    historical_summary_blob.upload_from_filename( HistSummaryFilePath )
    os.remove( HistSummaryFilePath )

    # upload historical data
    hist_prev_cloud_name = f"{group}-historical-prevalence.csv"
    HistPrevCloudPath = f"{CloudPathRoot}/group-{group}/{hist_prev_cloud_name}"
    historical_data_blob = bucket.blob( HistPrevCloudPath )
    historical_data_blob.upload_from_filename( HistPrevFilePath )

    ''' RUN SIMULATIONS '''

    # specified MDA coverage values
    for MDA_Cov in [ 0.6, 0.7, 0.8, 0.9 ]:

        try:

            # annual, biannual
            for mda_type, mda_vector in mda_vectors.items():

                mda_lists = [ mda_vector[i:j] for i, j in itertools.combinations( range( len( mda_vector ) + 1 ), 2 ) ]

                for mda_list in mda_lists:

                    # mkdir -p group dirs
                    output_file_path = f"group-{group}/coverage-{MDA_Cov}/mdatype-{mda_type}"
                    output_dir = f"{FilePathRoot}/output/{output_file_path}"
                    pathlib.Path( output_dir ).mkdir( parents = True, exist_ok = True )

                    # set up cloud storage path
                    cloud_dir = f"{CloudPathRoot}/{output_file_path}"

                    input_dir = f"{FilePathRoot}/input/{output_file_path}"
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

                    # output CSV file paths
                    prev_csv_file_name = f"{file_name_root}-prev.csv"
                    PrevFilePath = f"{output_dir}/{prev_csv_file_name}"
                    PrevCloudPath = f"{cloud_dir}/{prev_csv_file_name}"

                    infect_csv_file_name = f"{file_name_root}-infect.csv"
                    InfectFilePath = f"{output_dir}/{infect_csv_file_name}"
                    InfectCloudPath = f"{cloud_dir}/{infect_csv_file_name}"

                    print( f"=== Running:\n\tGroup: {group}\n\tMDA_Cov {MDA_Cov}\n\tMDA type {mda_type}\n\tinput {MDAFilePath}\n\toutput {PrevFilePath}\n" )

                    # run the simulation
                    Trachoma_Simulation(
                        BetFilePath=BetFilePath,
                        MDAFilePath=MDAFilePath,
                        PrevFilePath=PrevFilePath,
                        InfectFilePath=InfectFilePath,
                        SaveOutput=False,
                        OutSimFilePath=None,
                        InSimFilePath=InSimFilePath,
                        MDA_Cov=MDA_Cov
                    )

                    # upload prevalence file
                    prev_blob = bucket.blob( PrevCloudPath )
                    prev_blob.upload_from_filename( PrevFilePath )

                    # upload infection file
                    infect_blob = bucket.blob( InfectCloudPath )
                    infect_blob.upload_from_filename( InfectFilePath )

                    # remove the input file
                    os.remove( MDAFilePath )

                    ''' SUMMARISE OUTPUT DATA '''

                    # read in the simulation output
                    op_data = pd.read_csv( PrevFilePath )

                    # make a json file path to summarise it into
                    json_file_name = f"{file_name_root}-summary.json"
                    summary_json_path = f"{output_dir}/{json_file_name}"
                    summary_json_cloud_path = f"{cloud_dir}/{json_file_name}"

                    # summarise it in there
                    pd.DataFrame( {
                        'median': op_data.iloc[:, 2:].median(),
                        'lower': op_data.iloc[:, 2:].quantile(0.05),
                        'upper': op_data.iloc[:, 2:].quantile(0.95)
                    }).to_json( summary_json_path )

                    # upload summary file
                    summary_blob = bucket.blob( summary_json_cloud_path )
                    summary_blob.upload_from_filename( summary_json_path )

                    ''' REMOVE LOCAL OUTPUT FILES '''
                    os.remove( PrevFilePath )
                    os.remove( InfectFilePath )
                    os.remove( summary_json_path )

        except Exception as e:
            print( f"+++++ ERROR in group {group} MDA_Cov {MDA_Cov}: {str(e)}" )

    print( f"===== FINISHED RUNNING GROUP {group} =====" )

if len( groups ) > 1:
    print( f"***** FINISHED RUNNING GROUPS {groups} *****" )
