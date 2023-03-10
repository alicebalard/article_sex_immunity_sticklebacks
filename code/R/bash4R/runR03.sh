#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N runR
#$ -o /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/R/bash4R/run_R03.R.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/R/bash4R/run_R03.R.stderr

module load R/4.0.2

cd /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/R

Rscript /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/R/R03_prepObjectMethylkit_runInCLUSTER.R

