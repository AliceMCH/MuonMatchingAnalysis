#! /bin/bash

PID=$(ps ux | grep o2-dpl-run | grep session | head -n 1 | tr -s ' ' | cut -d' ' -f 2)

if [ x"$PID" != "x" ]; then
    echo "kill $1 $PID"
    kill $1 $PID
fi
