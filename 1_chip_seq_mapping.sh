#!/bin/bash

#####################################################
# ChIP-seq data: mapping and filtering
# (from fastq to filtered bam)
#####################################################

# mapping to hg19 or hg38

if [ $# -ne 3 ] ; then
   printf '\nUsage: %s genome input_fastq sample_name \n' "${0}";
   printf '       genome:      Specify genome - either hg19 or hg38\n';
   printf '       input_fastq:         Specify location of the input fastq file.\n';
   printf '       sample_name:         Provide sample name.\n';
   exit 1;
fi


module load bowtie/2.3.4
module load sambamba/0.5.5
module load htslib/1.8
module load samtools/1.8
module load bedtools/2.20.1 

genome="${1}"
input_fastq="${2}"
sample_name="${3}"

bowtie2 --very-sensitive --no-discordant -p 8 -x /data/groups/lab_bock/shared/resources/genomes/${genome}/indexed_bowtie2/${genome} --met-file ${sample_name}.aln_metrics.${genome}.txt ${input_fastq}  2> ${sample_name}.aln_rates.${genome}.txt | samtools view -S -b - | samtools sort -o ${sample_name}.trimmed.bowtie2.${genome}.bam -

sambamba markdup -t 8 -r --compression-level=0 ${sample_name}.trimmed.bowtie2.${genome}.bam ${sample_name}.trimmed.bowtie2.${genome}.filtered.nodups.nofilter.bam 2> ${sample_name}.dups_metrics.${genome}.txt
sambamba view -t 8 -f bam --valid -F "not unmapped and not (secondary_alignment or supplementary) and mapping_quality >= 30" ${sample_name}.trimmed.bowtie2.${genome}.filtered.nodups.nofilter.bam |sambamba sort -t 8 /dev/stdin -o ${sample_name}.trimmed.bowtie2.${genome}.filtered.bam

samtools index ${sample_name}.trimmed.bowtie2.${genome}.bam
samtools index ${sample_name}.trimmed.bowtie2.${genome}.filtered.bam

sambamba flagstat ${sample_name}.trimmed.bowtie2.${genome}.bam > ${sample_name}.trimmed.bowtie2.${genome}.flagstat
sambamba flagstat ${sample_name}.trimmed.bowtie2.${genome}.filtered.bam > ${sample_name}.trimmed.bowtie2.${genome}.filtered.flagstat


#bigWig
bedtools bamtobed -i ${sample_name}.trimmed.bowtie2.${genome}.filtered.bam | bedtools slop -i stdin -g /data/groups/lab_bock/shared/resources/genomes/${genome}/${genome}.chromSizes -s -l 0 -r 130 | /home/himrichova/1_Projects/0_scripts/fix_bedfile_genome_boundaries.py ${genome} | genomeCoverageBed -bg -g /data/groups/lab_bock/shared/resources/genomes/${genome}/${genome}.chromSizes -i stdin > ${sample_name}_${genome}.cov

awk 'NR==FNR{sum+= $4; next}{ $4 = ($4 / sum) * 1000000; print}' ${sample_name}_${genome}.cov ${sample_name}_${genome}.cov| sort -k1,1 -k2,2n > ${sample_name}_${genome}.normalized.cov

bedGraphToBigWig ${sample_name}_${genome}.normalized.cov /data/groups/lab_bock/shared/resources/genomes/${genome}/${genome}.chromSizes ${sample_name}_${genome}.bigWig

chmod 755 ${sample_name}_${genome}.bigWig
