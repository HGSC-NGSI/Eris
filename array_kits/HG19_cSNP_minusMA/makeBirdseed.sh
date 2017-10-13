#!/bin/bash

shopt -s nullglob

CONVERT_SCRIPT_DIR=`dirname $0`

perl ${CONVERT_SCRIPT_DIR}/csnpcsv2birdseed.pl \
     $1 ${CONVERT_SCRIPT_DIR}/probelist.txt

echo "$1 converted"
