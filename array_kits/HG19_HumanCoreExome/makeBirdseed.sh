#!/bin/bash

shopt -s nullglob

CONVERT_SCRIPT_DIR=`dirname $0`

perl ${CONVERT_SCRIPT_DIR}/csv2birdseed.pl \
     ${CONVERT_SCRIPT_DIR}/probelist.txt $1

echo "$1 converted"
