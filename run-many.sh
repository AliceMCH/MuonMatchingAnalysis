
#! /bin/bash

SCRIPT=$1
AO2DLIST=$2
OUTDIR="$3"

if [ -z "$OUTDIR" ]; then
    OUTDIR=AnalysisResults
fi

rm -f stop
rm -f ./AO2D.root
mkdir -p "$OUTDIR"


NFILES=$(cat "${AO2DLIST}" | wc -l)
I=1

while [ $I -le $NFILES ]; do

    if [ -e stop ]; then
        break;
    fi

    if [ -e "$OUTDIR"/AnalysisResults-${I}.root ]; then
        I=$((I+1))
        continue
    fi
    if [ -e "$OUTDIR"/done-${I} ]; then
        I=$((I+1))
        continue
    fi

    AO2DFILE=$(cat "${AO2DLIST}" | head -n $I | tail -n 1)
    echo "I=${I} => alien_cp \"alien://${AO2DFILE}\" file://./AO2D.root"
    rm -f ./AO2D.root && alien_cp "alien://${AO2DFILE}" file://./AO2D.root

    if [ -e stop ]; then
        break;
    fi

    echo "I=${I} => processing \"${AO2DFILE}\" ..."
    rm -f AnalysisResults.root FwdMatchMLCandidates.root MftDCA.root GlobalAlignment.root && $SCRIPT >& "$OUTDIR"/log-${I}.txt
    gzip -f "$OUTDIR"/log-${I}.txt
    echo "... done."

    if [ -e AnalysisResults.root ]; then
        touch "$OUTDIR"/done-${I}
        cp AnalysisResults.root "$OUTDIR"/AnalysisResults-${I}.root
    fi

    if [ -e FwdMatchMLCandidates.root ]; then
        cp FwdMatchMLCandidates.root "$OUTDIR"/FwdMatchMLCandidates-${I}.root
    fi

    if [ -e MftDCA.root ]; then
        cp MftDCA.root "$OUTDIR"/MftDCA-${I}.root
    fi

    if [ -e GlobalAlignment.root ]; then
        cp GlobalAlignment.root "$OUTDIR"/GlobalAlignment-${I}.root
    fi


    I=$((I+1))

    #break;

    #rm -f ./AO2D.root

done

echo ""
echo "Processing finished"
