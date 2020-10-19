'''
*Creating the iu-level.csv for Trachoma*

1) Work out the mean prevalence for all years (take the last month of the year 12-2015) for all groups based on the group files OutputPrev_groupXX

2) Using the IUcodeTrachomaGroupsEthiopia.csv for Ethiopia or if we're doing all areas IUcodeTrachomaGroups.csv work out all IUIDs (IUcodes) that each group encompasses

We then need to assemble the a file in the format as per example:

public/diseases/trachoma/iu-level.csv

 IUID = IUcode

Endemicity = Unknown
'''

import re
import sys
import glob
import json
import pandas as pd

# TODO paramaterize paths
# TODO parameterize groups_to_use

groupToIus = {}
iuGroupMappings = json.load( open("../../ntd-simulator/src/pages/components/simulator/helpers/iuGroupMapping.json"))
for k, v in iuGroupMappings.items():
    gid = str(v)
    if not gid in groupToIus:
        groupToIus[ gid ] = []
    if not k in groupToIus[ gid ]:
        groupToIus[ gid ].append( k )

output_data = []

output_columns = [ 'IUID', 'Endemicity' ]
output_columns.extend( [ f"Prev_Year{x}" for x in range( 2000, 2020 ) ] )

output_prev_group_file_path = 'data/merged'
pattern = "^" + output_prev_group_file_path + "/OutputPrev_group(?P<gid>[0-9]{1,3}).csv$"

groups_to_use = [ 99, 103, 124, 126, 128, 132, 142 ]

for PrevFilePath in glob.glob(f"{output_prev_group_file_path}/OutputPrev_group[0-9]*.csv"):

    gid = re.match( pattern, PrevFilePath ).group('gid')

    if not int( gid ) in groups_to_use:
        continue

    sys.stderr.write( f"processing file {PrevFilePath} for group {gid}\n" )

    data = pd.read_csv( PrevFilePath )
    cols_we_want = [ c for c in list( data ) if c[0:2] == '12' ]
    medians = data[ cols_we_want ].median().to_dict()

    ius = groupToIus[ gid ]
    for iu in groupToIus[ gid ]:
        sys.stderr.write( f"\tgenerating row for iu {iu} in group {gid}\n" )
        new_row = [ iu, 'Unknown' ]
        for column in output_columns[ 2: ]:
            year = column[9:]
            median_key = f"12-{year}"
            if median_key in medians:
                new_row.append( medians[ median_key ] )
            else:
                new_row.append( 'NA' )

        output_data.append( new_row )

df = pd.DataFrame( output_data, columns = output_columns )
print( df.to_csv( index = None ) )
