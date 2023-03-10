#!/bin/bash                                                                                                                                                         
#$ -pe smp 15
#$ -l h_vmem=5G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y                                                                                                                                             
#$ -N runR05
#$ -o /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/R/bash4R/run_R05.R.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/R/bash4R/run_R05.R.stderr

module load R/4.0.2
cd /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/R
Rscript /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/R/R05_getDifferentialMethylation_runInCLUSTER.R
