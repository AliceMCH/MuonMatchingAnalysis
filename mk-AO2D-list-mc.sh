#! /bin/bash

if [ $# -lt 2 ]; then
    echo "usage: mk-AO2D-list-mc.sh year prod pass"
    exit 1
fi

YEAR=$1
PROD=$2
PASS=$3

#echo "alien_find /alice/sim/$YEAR/$PROD/$PASS \"*/AOD/*/AO2D.root\""
#alien_find /alice/sim/$YEAR/$PROD/$PASS "*/AOD/*/AO2D.root"

#echo "alien_find /alice/sim/$YEAR/$PROD/$PASS \"*/AOD/*/AO2D.root\""
alien_find /alice/sim/$YEAR/$PROD/$PASS "*/*/AO2D.root"
