#!/bin/bash -e
JSON_FILE=$(realpath "$1")
RESULTS=$(realpath "$2")
OUTPUT=$(realpath "$3")

DIR="$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )"
INSTALL="$( dirname "$DIR" )"

cd $RESULTS

rootline="${DIR}/macros/runDraw.C(\"$JSON_FILE\",\"${OUTPUT}\",\"${INSTALL}\")"
#echo root -l -b -q $rootline
root -l -b -q $rootline
