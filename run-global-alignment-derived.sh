#! /bin/bash

OPTION="-b --configuration json://config.json --pipeline muon-global-alignment:1 -b"

AO2DFILE=$1
if [ -z "$AO2DFILE" ]; then
    AO2DFILE=AO2D.root
fi

echo "o2-analysis-dq-muon-global-alignment ${OPTION} --aod-file \"$AO2DFILE\""
o2-analysis-dq-muon-global-alignment ${OPTION} --aod-file "$AO2DFILE"
