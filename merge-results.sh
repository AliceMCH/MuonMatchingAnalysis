
#! /bin/bash

DIR="$1"
if [ -z "$DIR" ]; then
    DIR="AnalysisResults"
fi

WD=$(dirname $0)

$WD/merge-root-files.sh $DIR "AnalysisResults"
$WD/merge-root-files.sh $DIR "FwdMatchMLCandidates"
$WD/merge-root-files.sh $DIR "MftDCA"

echo ""
echo "Merging done"
