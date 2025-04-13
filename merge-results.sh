
#! /bin/bash

DIR="$1"
if [ -z "$DIR" ]; then
    DIR="AnalysisResults"
fi

hadd -f ${DIR}/AnalysisResultsFull.root $DIR/AnalysisResults-*.root
