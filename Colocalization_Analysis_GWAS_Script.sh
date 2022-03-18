#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=20GB
#PBS -l walltime=04:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

#==============================
# Colocalization_Analysis script
#==============================

CurrDir=`pwd`

CodeExec=$CurrDir'/Colocalization_Analysis_GWAS_Script.R'
ConfigFile=$CurrDir'/configfile.yaml'

Rscript $CodeExec $ConfigFile


