#! /bin/bash

if [ $# -lt 1 ]; then
    echo "usage: mk-AO2D-list-derived.sh joblist"
    exit 1
fi

JOBLIST=$1

I=0
while IFS= read -r -d',' job
do
    
    echo "alien://${job}/AO2D.root"

    I=$((I+1))

done < ${JOBLIST}
