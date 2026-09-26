#! /bin/bash

OPTION="-b --configuration json://configuration.json"

o2-analysis-track-to-collision-associator ${OPTION} | \
    o2-analysis-trackselection ${OPTION} | \
    o2-analysis-track-propagation ${OPTION} | \
    o2-analysis-mm-track-propagation ${OPTION} | \
    o2-analysis-fwdtrackextension ${OPTION} | \
    o2-analysis-fwdtrack-to-collision-associator ${OPTION} | \
    o2-analysis-multcenttable ${OPTION} | \
    o2-analysis-event-selection-service ${OPTION} | \
    o2-analysis-dq-mft-mch-matcher ${OPTION} --aod-file AO2D.root --aod-writer-json writer.json --shm-segment-size 75000000000
