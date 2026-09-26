#! /bin/bash

if [ $# -lt 2 ]; then
    echo "usage: mk-AO2D-list-mc.sh year prod pass"
    exit 1
fi

YEAR=$1
PROD=$2
PASS=$3

alien_find /alice/data/$YEAR/$PROD "*/$PASS/*/*/*/AO2D.root"
