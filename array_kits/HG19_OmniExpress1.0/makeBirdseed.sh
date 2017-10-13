#!/bin/bash

shopt -s nullglob

CONVERT_SCRIPT_DIR=`dirname $0`

python ${CONVERT_SCRIPT_DIR}/makeBirdseed.py \
       ${CONVERT_SCRIPT_DIR}/probelist.txt $1

echo "$1 converted"
