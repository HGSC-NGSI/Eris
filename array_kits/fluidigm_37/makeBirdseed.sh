#!/bin/bash

inFile=$1
echo "Directory: $PWD"
echo "Raw array: $inFile"

CONVERT_SCRIPT_DIR=`dirname $0`

python3 ${CONVERT_SCRIPT_DIR}/csv2birdseed.py \
    ${CONVERT_SCRIPT_DIR}/probelist.txt \
    "$inFile" \
    ${CONVERT_SCRIPT_DIR}/RSID_RUID_lookup.txt

echo "converted $(basename ${inFile})"
echo
