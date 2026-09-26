#! /bin/bash

if [ $# -lt 2 ]; then
    echo "usage: get-AO2D-mc.sh year prod pass index"
    exit 1
fi

YEAR=$1
PROD=$2
PASS=$3
INDEX=$4

./mk-AO2D-list-mc.sh $YEAR $PROD $PASS > _temp.txt

LINE=$(cat _temp.txt | sed -n ${INDEX}p)

echo "AO2D: $LINE"

if [ -z "$LINE" ]; then exit 1; fi

rm -f ./AO2D.root
alien_ls -l "alien://$LINE"
alien_cp "alien://$LINE" file://./AO2D.root
