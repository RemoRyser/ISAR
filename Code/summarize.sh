#!/bin/bash


#$ -S /bin/bash

#$ -M remo.ryser@idiv.de
#$ -m a
#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.out
#$ -e /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.err
#$ -cwd

#$ -l h_rt=2:00:00
#$ -l h_vmem=1G

#$ -binding linear:1





INPUTDIR=$1
WEB=$2
NAME=$3

cd "$INPUTDIR" || exit

OUTFILEGLOBAL=SummaryGlobal.out
OUTFILEMASS=SummaryMass.out

mkdir "/data/idiv_brose/remo/Patchdistance/Simulation/output/$NAME$WEB"

echo 'landscape, rep, number.of.spp, number.of.plants, number.of.consumers, number.of.patch, rng.seed, max.emigr.rate.plants, shape.emigr.rate.plants, max.emigr.rate.consumers, shape.emigr.rate.consumers, D_0, theta, eps, mean.patch.dist, sd.patch.dist, mean.nn.dist, sd.nn.dist, mean.con.rgg, sd.con.rgg, mean.con.rgg.plants, sd.con.rgg.plants, mean.con.rgg.consumers, sd.con.rgg.consumers, Ricker' > $OUTFILEGLOBAL
for file in global_*.out ; do

N=$(echo "$file" | grep -Eo '[0-9]+')


ANALYSIS=NULL,
REP=$N,
dat=$(sed '1d' global_"${N}".out)

echo "$ANALYSIS $REP $dat"  >> "$OUTFILEGLOBAL"

done





echo 'landscape, rep, patch, species, body.mass, int.biomass, mean.biomass,  biomass.variance, tot.mean.biomass, tot.biomass.variance, if.basal.spp, biomass.tend-20k, biomass.tend-10k,biomass.tend, spatial.connectance, Mean.net.growth1, Mean.net.growth2, net.growth.sd, Area'  > $OUTFILEMASS


for file in mass_*.out ; do

N=$(echo "$file" | grep -Eo '[0-9]+')


sed '1d' mass_"${N}".out > temp.out

ANALYSIS=NULL,
REP=$N,

while  read -r p; do


echo "$ANALYSIS $REP $p" >> $OUTFILEMASS

done < temp.out

rm temp.out

cp web_*.out "/data/idiv_brose/remo/Patchdistance/Simulation/output/$NAME$WEB/"
cp params_*.out "/data/idiv_brose/remo/Patchdistance/Simulation/output/$NAME$WEB/"


done

cp  SummaryGlobal.out "/data/idiv_brose/remo/Patchdistance/Simulation/output/$NAME$WEB/SummaryGlobal_$WEB.out"
cp  SummaryMass.out "/data/idiv_brose/remo/Patchdistance/Simulation/output/"$NAME$WEB"/SummaryMass_$WEB.out"

