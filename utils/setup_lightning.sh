#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <simulation name>"
    exit 1
fi

simname="$1"
scratchpath="${HOME}/scratch/kilonovae"

codepath="${HOME}/software/thunderstorm"
lcodepath="${codepath}/lightning"

simpath="${scratchpath}/${simname}"
lsimpath="${simpath}/lightning"

mkdir -p "${simpath}"
mkdir -p "${simpath}/resources"
mkdir -p "${lsimpath}"

mkdir -p "${lsimpath}/build"
mkdir -p "${lsimpath}/config"
mkdir -p "${lsimpath}/run"
mkdir -p "${lsimpath}/data"

cp "${codepath}/LICENSE" "${simpath}/"
cp "${codepath}/resources/*" "${simpath}/resources/"
cp "${lcodepath}/build/lightning" "${lsimpath}/build/"
cp "${lcodepath}/config/config.ini" "${lsimpath}/config/"
cp "${lcodepath}/run/run.sh" "${lsimpath}/run/"

