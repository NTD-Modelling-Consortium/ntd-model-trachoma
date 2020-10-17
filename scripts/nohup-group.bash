#!/usr/bin/env bash

nohup pipenv run python scripts/run_group.py $1 > run_group.out 2> run_group.err < /dev/null &
