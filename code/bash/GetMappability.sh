#!/bin/bash
# Bash script to output number of reads total and mapped, with and without sex chromosome + 
unknown chromosome
# NB: Bash can't deal with non integer, the percentage mappability can be calculated a posteriori e.g. in R
#$ -N mappability
#$ -o /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/bash/run_mappability.stdout
#$ -e /data/SBCS-EizaguirreLab/Alice/SexImmunity/article_sex_immunity_sticklebacks/code/bash/run_mappability.stderr
#$ -cwd
#$ -V
#$ -l h_rt=240:00:00

module load samtools/1.9

OUTDIR=/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/04BSBolt_methCall/BSBolt/Alignments
cd $OUTDIR

# Create file:
echo "filename" "Total_reads" "Mapped_read" "Total_reads_minusSexAndUnknownChr" "Mapped_reads_minusSexAndUnknownChr" > Report_mapping_efficiency.txt

for file in `ls $OUTDIR/*fastq.gz.bam.sorted.bam`
do
    # Number of mapped reads excluding those marked as secondary (some other probable hits for the read) or supplementary (i.e., multimappers or chimeric entries)
    MAPPED=$(samtools view -c -F 0x904 $file)

    # Count total number of reads by read name (avoid secondary):
    TOTAL=$(samtools view -c -F 2048 $file)

    # Number of mapped reads WITHOUT sex & unknown chromosomes:
    MAPPED2=$(samtools idxstats $file | cut -f 1 | grep -vE 'Gy_chrUn|Gy_chrXIX' | xargs samtools view -b $file -c -F 0x904)

    # Count total number of reads WITHOUT sex & unknown chromosomes:
    TOTAL2=$(samtools idxstats $file | cut -f 1 | grep -vE 'Gy_chrUn|Gy_chrXIX' | xargs samtools view -b $file -c -F 2048)

    # output in one line in the file
    echo $file $TOTAL $MAPPED $TOTAL2 $MAPPED2 >> Report_mapping_efficiency.txt
done
