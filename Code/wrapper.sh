#!/bin/bash

if [[ -z $2 ]] ; then
echo "usage: bash $0 jobname number.of.landscapes number.of.foodwebs"
exit 1
fi


LANDSCAPE=$2


INPUTDIR=$(dirname $0)

for WEB in $(seq 1 40)

do


NAME=$1


JOB_ID=$(
qsub \
-terse \
${var:+-v VAR=$var} \
-t 1-100 \
-N $NAME \
$(dirname $0)/submit.sh $NAME $WEB $INPUTDIR $LANDSCAPE| cut -d. -f1)



qsub -hold_jid $JOB_ID $(dirname $0)/summarize.sh /work/$USER/$NAME/$WEB $WEB $NAME


done



