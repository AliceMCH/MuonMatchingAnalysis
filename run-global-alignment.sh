#! /bin/bash

OPTION="-b --configuration json://config.json --pipeline track-propagation:1,muon-global-alignment:1 -b"
#QA_PROF="--child-driver 'valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all  --track-origins=yes --verbose --log-file=valgrind-out.txt'"
QA_PROF=""


echo "o2-analysis-dq-muon-global-alignment ${OPTION} ${QA_PROF} | o2-analysis-event-selection ${OPTION} | o2-analysis-fwdtrackextension ${OPTION} | o2-analysis-track-propagation ${OPTION} | tracksExtraConverter-001-to-002 ${OPTION} | o2-analysis-timestamp ${OPTION}"

o2-analysis-dq-muon-global-alignment ${OPTION} ${QA_PROF} | o2-analysis-event-selection ${OPTION} | o2-analysis-fwdtrackextension ${OPTION} | o2-analysis-track-propagation ${OPTION} | o2-analysis-timestamp ${OPTION} --aod-file AO2D.root --aod-writer-json writer-global-alignment.json
#o2-analysis-dq-qa-matching-dev ${OPTION} ${QA_PROF} | o2-analysis-event-selection ${OPTION} | o2-analysis-fwdtrackextension ${OPTION} | o2-analysis-track-propagation ${OPTION} | o2-analysis-timestamp ${OPTION} --aod-file AO2D.root --aod-writer-json writer-global-alignment.json

#o2-analysis-dq-muon-global-alignment ${OPTION} ${QA_PROF} | o2-analysis-event-selection ${OPTION} | o2-analysis-fwdtrackextension ${OPTION} | o2-analysis-track-propagation ${OPTION} | o2-analysis-tracks-extra-v002-converter ${OPTION} | o2-analysis-timestamp ${OPTION}

#o2-analysis-dq-muon-global-alignment ${OPTION} ${QA_PROF} --aod-writer-json writer-global-alignment.json | o2-analysis-event-selection ${OPTION} | o2-analysis-fwdtrackextension ${OPTION} | o2-analysis-track-propagation ${OPTION} | o2-analysis-timestamp ${OPTION} --shm-segment-size 75000000000
