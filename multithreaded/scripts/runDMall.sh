#!/bin/bash -e

INPUT=$(realpath "$1")
JSON_FILE=$(basename "$INPUT")

DIR="$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )"
INSTALL="$( dirname "$DIR" )"

mkdir -p results
if [ -f .rootrc ]; then
    cp .rootrc results
fi

cd results
for sign in muplus muminus
do
	cp $INPUT $JSON_FILE.$sign
	# Define the key to be replaced and the new value
	KEY="sign"
	NEW_VALUE=$sign
	# Use sed to replace the old value with the new value
	sed -i "s/\"$KEY\":.*\".*\"/\"$KEY\": \"$NEW_VALUE\"/" "$JSON_FILE.$sign"

	rootline="${DIR}/macros/run_D_M_Every.C(\"$JSON_FILE.$sign\",\"${INSTALL}\")"
	echo root -l -b -q $rootline
	root -l -b -q $rootline
done
