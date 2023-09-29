#!/bin/bash

if [[ -z $2 ]] ; then
echo "usage: bash $0 jobname number.of.foodwebs"
exit 1
fi




NAME=$1
Iterations=1

for WEB in $(seq 1 $2)

do

LIST=""

SEED=$WEB

JOB_ID=$(
qsub \
-terse \
${var:+-v VAR=$var} \
-t 1-$Iterations \
-N $NAME$WEB \
$(dirname $0)/submit.sh $SEED $NAME| cut -d. -f1)



LIST+="$JOB_ID,"


done

qsub -hold_jid ${LIST%%,} $(dirname $0)/summarize.sh /work/$USER $WEB $NAME



