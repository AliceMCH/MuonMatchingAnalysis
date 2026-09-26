#! /bin/bash

OPTION="-b --configuration json://config-matching-qa.json --pipeline track-propagation:1,qa-matching:1 -b"
o2-analysis-track-to-collision-associator ${OPTION} | \
    o2-analysis-trackselection ${OPTION} | \
    o2-analysis-track-propagation ${OPTION} | \
    o2-analysis-mm-track-propagation ${OPTION} | \
    o2-analysis-fwdtrack-to-collision-associator ${OPTION} | \
    #o2-analysis-dq-global-muon-matcher ${OPTION} | \
    o2-analysis-fwdtrackextension ${OPTION} | \
    o2-analysis-multcenttable ${OPTION} | \
    o2-analysis-event-selection-service ${OPTION} | \
    #o2-analysis-dq-mft-mch-matcher ${OPTION} --aod-writer-json writer.json | \
    o2-analysis-dq-qa-matching ${OPTION} --shm-segment-size 75000000000 --aod-file AO2D.root # --aod-writer-json writer-matching-qa.json
