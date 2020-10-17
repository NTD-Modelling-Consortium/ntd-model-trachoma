#!/usr/bin/env bash

nohup pipenv run python scripts/run_group.py "$1" > "run_group_$1.out" 2> "run_group_$1.err" < /dev/null &
echo "--- started background group run for group $1 on host $( hostname )"
