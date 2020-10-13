#!/usr/bin/env bash

DIR=data/merged

for group in $( ls -1 ${DIR}/InputBet_group*.csv | sed 's/.*_group\([0-9]*\)\.csv/\1/g' ) ; do
    [[ -f "${DIR}/OutputVals_group${group}.p" ]] && python run_group.py ${group} || echo "OutputVals file not found for group $group"
done
