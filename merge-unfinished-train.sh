#! /bin/bash

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
    rootfiles=$(alien_find ${job} "/*/AnalysisResults.root" | grep -v -e Stage)

    fileindex=1
    for f in ${rootfiles}; do
        echo "  $f"
        alien_cp alien://"$f" file://./${jobdir}/AnalysisResults-${fileindex}.root
        fileindex=$((fileindex+1))
    done
    ./merge-root-files.sh ${jobdir} AnalysisResults
    cp ${jobdir}/AnalysisResultsFull.root _tmp/AnalysisResults-${jobindex}.root

    jobindex=$((jobindex+1))
done < ${joblist}
