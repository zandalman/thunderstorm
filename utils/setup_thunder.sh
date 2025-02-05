#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <simulation name> <pp name>"
    exit 1
fi

simname="$1"
ppname="$2"
scratchpath="${HOME}/scratch/kilonovae"

codepath="${HOME}/software/thunderstorm"
tcodepath="${codepath}/thunder"

simpath="${scratchpath}/${simname}"
tsimpath="${simpath}/${ppname}"

mkdir -p "${simpath}"
mkdir -p "${tsimpath}"

mkdir -p "${tsimpath}/build"
mkdir -p "${tsimpath}/config"
mkdir -p "${tsimpath}/run"
mkdir -p "${tsimpath}/data"
mkdir -p "${tsimpath}/data/hist"

cp "${tcodepath}/build/thunder" "${tsimpath}/build/"
cp "${tcodepath}/config/config.ini" "${tsimpath}/config/"
cp "${tcodepath}/run/run.sh" "${tsimpath}/run/"

