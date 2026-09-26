
#! /bin/bash

rm -f stop
rm -f ./AO2D.root
#rm -rf AnalysisResults
mkdir -p AnalysisResults

#AO2DLIST=AO2D_list-LHC24l7.txt
#AO2DLIST=AO2D_list-LHC22p.txt
#AO2DLIST=AO2D_list-LHC23h-apass4_skimmed.txt
#AO2DLIST=AO2D_list-LHC23zk.txt
#AO2DLIST=AO2D_list-LHC24am.txt
#AO2DLIST=AO2D_list-LHC24aq-apass1_muon_matching.txt

# 2024 ppRef with orfficial alignment
#AO2DLIST=AO2D_list-LHC24aq-apass1.txt
#AO2DLIST=AO2D_list-LHC24aq-apass1-559408-0000.txt

# 2024 ppRef with new MFT geometry
AO2DLIST=AO2D_list-LHC24aq-apass1_muon_matching3.txt
#AO2DLIST=AO2D_list-LHC24aq-apass1_muon_matching3-559408-0000.txt

# 2024 B-off with new MFT geometry
#AO2DLIST=AO2D_list-LHC24ad-apass4_muon_matching.txt

# 2025 OO apass2
#AO2DLIST=AO2D_list-LHC25ae-apass2.txt

#AO2DLIST=AO2D_list-LHC25c3b.txt

#AO2DLIST=AO2D_list-LHC25i4.txt

#AO2DLIST=AO2D_list-LHC24ar_apass1.txt
#AO2DLIST=AO2D_list-LHC24ar_apass3.txt



NFILES=$(cat "${AO2DLIST}" | wc -l)
I=1

while [ $I -le $NFILES ]; do

    if [ -e stop ]; then
        break;
    fi

    if [ -e AnalysisResults/AnalysisResults-${I}.root ]; then
        I=$((I+1))
        continue
    fi
    if [ -e AnalysisResults/done-${I} ]; then
        I=$((I+1))
        continue
    fi

    AO2DFILE=$(cat "${AO2DLIST}" | head -n $I | tail -n 1)
    echo "I=${I} => alien_cp \"alien://${AO2DFILE}\" file://./AO2D.root"
    rm -f ./AO2D.root && alien_cp "alien://${AO2DFILE}" file://./AO2D.root

    echo "I=${I} => processing \"${AO2DFILE}\" ..."
    rm -f AnalysisResults.root && ./run-qa.sh >& AnalysisResults/log-${I}.txt
    gzip -f AnalysisResults/log-${I}.txt
    echo "... done."

    if [ -e AnalysisResults.root ]; then
        touch AnalysisResults/done-${I}
        cp AnalysisResults.root AnalysisResults/AnalysisResults-${I}.root
    fi

    I=$((I+1))

    #break;

    #rm -f ./AO2D.root

done
