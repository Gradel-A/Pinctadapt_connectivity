#!/usr/bin/bash

#declare variables used in the script
DATADIRECTORY=/analyses_article1
DATAINPUT=$DATADIRECTORY/02_data
DATAOUTPUT=$SCRATCH/analyses_article1/admixture
SCRIPT=$DATADIRECTORY/00_scripts/admixture
HEADER=$DATADIRECTORY/00_scripts/headerP.txt
NCPU=28
ADMIXTURENV=". /admixture/1.3.0/env.sh"


#create the output files
mkdir -p $SCRIPT
mkdir -p $DATAOUTPUT

#go to the folder containing the data in your project folder
cd $SCRATCH/06_freebayes/partial_filter
ls -d $DATADIRECTORY/02_data/*l.bed> $DATADIRECTORY/01_infofiles/input_files_admixture.txt

NAME='cat /analyses_article1/01_infofiles/input_files_admixture.txt'


for KANC in $(seq 1 10)
do
for FILE in $($NAME) #list all the files with the fna extension and store it in the variable FILE
do
        cp $HEADER $SCRIPT/run_${KANC}_${FILE##*/}.sh ; # copy the header file in a new script file called run_{genome_part}.sh
        echo "#PBS -N run_${KANC}_${FILE##*/}" >> $SCRIPT/run_${KANC}_${FILE##*/}.sh ;
        echo "#PBS -o $DATADIRECTORY/98_logfiles/10_run_${KANC}_${FILE##*/}.log" >> $SCRIPT/run_${KANC}_${FILE##*/}.sh ;
        echo "cd $DATAOUTPUT" >> $SCRIPT/run_${KANC}_${FILE##*/}.sh ; # append the echoed line in the script file (we go to the data folder of our project)
        echo "$ADMIXTURENV" >> $SCRIPT/run_${KANC}_${FILE##*/}.sh ; # preparing the run
        echo "admixture --cv --seed=25 ${FILE} ${KANC} -j28 | tee log_admixture_${KANC}_${FILE##*/}.log" >> $SCRIPT/run_${KANC}_${FILE##*/}.sh ;
        qsub $SCRIPT/run_${KANC}_${FILE##*/}.sh ; # append the echoed line in the script file (here we ask to submit our script to the PBS claculation nodes for execution)
done ; # we finish the loops
done ;