import re
import glob
import pathlib
import pandas as pd

# TODO parameterize paths in *argv

for PrevFilePath in glob.glob('data/merged/OutputPrev_group[0-9]*.csv'):

    output_prevalence = pd.read_csv( PrevFilePath )

    gid = re.match( "^data/merged/OutputPrev_group(?P<gid>[0-9]{2,3}).csv$", PrevFilePath ).group('gid')

    # FIXME
    output_dir = f"../../ntd-simulator/public/diseases/trachoma/data/group-{gid}"
    pathlib.Path( output_dir ).mkdir( parents = True, exist_ok = True )

    summary = pd.DataFrame( {
        'median': output_prevalence.iloc[:, 2:].median(),
        'lower': output_prevalence.iloc[:, 2:].quantile(0.05),
        'upper': output_prevalence.iloc[:, 2:].quantile(0.95)
    })

    output_file = f"{output_dir}/{gid}-historical-prevalence-summary.json"
    summary.to_json( output_file )
