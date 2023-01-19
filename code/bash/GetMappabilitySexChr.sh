#!/bin/bash
# Bash script to output number of reads total and mapped on the sex chromosome + unknown chromosome
#$ -N mappability2
#$ -o /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/bash/run_mappability2.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/bash/run_mappability2.stderr
#$ -cwd
#$ -V
#$ -l h_rt=1:00:00

module load samtools/1.9

OUTDIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/Alignments
cd $OUTDIR

# Create file:
echo "filename" "Total_reads" "Mapped_read_sexchr" "Mapped_reads_unknownchr" > Report_mapping_efficiencySexUn.txt

for file in `ls $OUTDIR/*fastq.gz.bam.sorted.bam`
do
    # Count total number of reads by read name (avoid secondary):
    TOTAL=$(samtools view -c -F 2048 $file)
    
    # Number of mapped reads on sex chromosome:
    MAPPEDSEX=$(echo 'Gy_chrXIX' | xargs samtools view -b $file -c -F 0x904)

    # Number of mapped reads on unknown chromosome:
    MAPPEDUN=$(echo 'Gy_chrUn' | xargs samtools view -b $file -c -F 0x904)

    # output in one line in the file
    echo $file $TOTAL $MAPPEDSEX $MAPPEDUN >> Report_mapping_efficiencySexUn.txt
done
