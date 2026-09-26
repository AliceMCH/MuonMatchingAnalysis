
#! /bin/bash

DIR="$1"
PREFIX="$2"

FLIST=$(ls -1 $DIR/${PREFIX}-?.root $DIR/${PREFIX}-??*.root)

if [ -z "$FLIST" ]; then
    exit 1
fi

if [ -e ${DIR}/${PREFIX}Full.root ]; then
    cp ${DIR}/${PREFIX}Full.root ${DIR}/${PREFIX}Full.root.bak
fi

hadd -f ${DIR}/${PREFIX}Full.root $FLIST
