#!/bin/bash



if [[ -z $1 ]] ; then
echo "usage: bash $0 replicates"
exit 1
fi




OUTPUTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/Newoutput/"

for REP in $(seq 1 $1)

do



SEED=$REP


export SEED
export OUTPUTDIR


./test SEED OUTPUTDIR



done

