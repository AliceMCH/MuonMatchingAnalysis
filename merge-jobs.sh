#! /bin/bash

export SCRIPTDIR=$(readlink -f $(dirname $0))

#/alice/cern.ch/user/a/alihyperloop/jobs/0562/hy_5621583/0021

joblist=$1

#jobs=$(alien_ls /alice/cern.ch/user/a/alihyperloop/jobs/${train})

#echo "$jobs"

rm -rf _tmp
jobindex=1
while IFS= read -r -d',' job
do
    jobdir=_tmp/${jobindex}
    rm -rf ${jobdir}
    mkdir -p ${jobdir}
    echo $job
    alien_cp alien://${job}/AnalysisResults.root file://./_tmp/AnalysisResults-${jobindex}.root

    jobindex=$((jobindex+1))
done < ${joblist}

$SCRIPTDIR/merge-root-files.sh _tmp AnalysisResults
