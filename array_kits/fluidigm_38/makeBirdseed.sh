#!/bin/bash

inFile=$1
echo "Directory: $PWD"
echo "Raw array: $inFile"

CONVERT_SCRIPT_DIR=`dirname $0`

python ${CONVERT_SCRIPT_DIR}/csv2birdseed.py \
    ${CONVERT_SCRIPT_DIR}/probelist.txt \
    "$inFile" \
    ${CONVERT_SCRIPT_DIR}/RSID_RUID_lookup.txt

echo "Converted $(basename $inFile)"
echo
