#! /bin/bash

echo "ALIBUILD_WORK_DIR: ${ALIBUILD_WORK_DIR}"

#(cd ${ALIBUILD_WORK_DIR}/BUILD/O2Physics-latest/O2Physics && ninja stage/bin/o2-analysis-dq-task-muon-dca && ninja stage/bin/o2-analysis-dq-table-maker)

#(cd ${ALIBUILD_WORK_DIR}/BUILD/O2Physics-latest/O2Physics && ninja stage/bin/o2-analysis-dq-mft-mch-matcher)
#exit

#(cd ${ALIBUILD_WORK_DIR}/BUILD/O2Physics-latest/O2Physics && ninja PWGDQ/install) && \
(cd ${ALIBUILD_WORK_DIR}/BUILD/O2Physics-latest/O2Physics && ninja PWGDQ/Tasks/install) || exit 1
(cd ${ALIBUILD_WORK_DIR}/BUILD/O2Physics-latest/O2Physics && ninja Common/Tasks/install) || exit 1
(cd ${ALIBUILD_WORK_DIR}/BUILD/O2Physics-latest/O2Physics && ninja Common/TableProducer/install) || exit 1
(cd ${ALIBUILD_WORK_DIR}/BUILD/O2Physics-latest/O2Physics && ninja PWGDQ/TableProducer/install) || exit 1
