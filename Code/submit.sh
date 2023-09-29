#!/bin/bash
 
#$ -S /bin/bash

##$ -M remo.ryser@idiv.de
#$ -m a
#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.out
#$ -e /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.err
#$ -cwd

#$ -l h_rt=2:00:00
#$ -l h_vmem=4G

#$ -binding linear:1



WEB=$2
LANDSCAPE=NULL
SEED=$SGE_TASK_ID
NAME=$1

module load gsl/1.16-3
module load sundials/2.7.0-1

OUTPUTDIR="/work/$USER/$NAME/$WEB/"
INPUTDIR="$3"
mkdir -p "$OUTPUTDIR"


export OUTPUTDIR
export NAME
export LANDSCAPE
export SEED
export WEB
export INPUTDIR

./simulation




