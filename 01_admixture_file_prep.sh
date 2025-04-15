#!/usr/bin/env bash
#PBS -N admixture_art
#PBS -q mpi
#PBS -l ncpus=28
#PBS -l mem=60gb
#PBS -l walltime=48:00:00
#PBS -o /admixture_fileprep.log

#Global variables
DATADIRECTORY=/analyses_article1
DATAINPUT=$DATADIRECTORY/02_data
DATAINFO=$DATADIRECTORY/01_infofiles
DATAOUTPUT=$SCRATCH/analyses_article1/admixture

PLINK=". /plink/1.9/env.sh"
NCPU=28

mkdir -p $DATAOUTPUT

# first we will transforme the vcf file into the bed format require by admixture, select the good populations and filtrate again
cd $DATAINPUT

$PLINK

plink --allow-extra-chr --vcf $DATAINPUT/filtered_data_hwe.vcf --make-bed --out $DATAINPUT/filtered_data_hwe

plink --allow-extra-chr --bfile $DATAINPUT/filtered_data_hwe --extract $DATAINFO/outliers.ID.txt --make-bed --out $DATAINPUT/outlier_dataset_all
plink --allow-extra-chr --bfile $DATAINPUT/filtered_data_hwe --exclude $DATAINFO/outliers.ID.txt --make-bed --out $DATAINPUT/neutral_dataset_all
plink --allow-extra-chr --bfile $DATAINPUT/filtered_data_hwe --keep $DATAINFO/indiv.pol.txt --extract $DATAINFO/outliers.ID.pol.txt --make-bed --out $DATAINPUT/outlier_dataset_pol
plink --allow-extra-chr --bfile $DATAINPUT/filtered_data_hwe --keep $DATAINFO/indiv.pol.txt --exclude $DATAINFO/outliers.ID.pol.txt --make-bed --out $DATAINPUT/neutral_dataset_pol

rm $SCRATCH/analyses_article1/tmp*

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' $DATAINPUT/outlier_dataset_all.bim > $DATAINPUT/outlier_dataset_all.bim.tmp
awk '{$1="0";print $0}' $DATAINPUT/neutral_dataset_all.bim > $DATAINPUT/neutral_dataset_all.bim.tmp
awk '{$1="0";print $0}' $DATAINPUT/outlier_dataset_pol.bim > $DATAINPUT/outlier_dataset_pol.bim.tmp
awk '{$1="0";print $0}' $DATAINPUT/neutral_dataset_pol.bim > $DATAINPUT/neutral_dataset_pol.bim.tmp

mv $DATAINPUT/outlier_dataset_all.bim.tmp $DATAINPUT/outlier_dataset_all.bim
mv $DATAINPUT/neutral_dataset_all.bim.tmp $DATAINPUT/neutral_dataset_all.bim
mv $DATAINPUT/outlier_dataset_pol.bim.tmp $DATAINPUT/outlier_dataset_pol.bim
mv $DATAINPUT/neutral_dataset_pol.bim.tmp $DATAINPUT/neutral_dataset_pol.bim

