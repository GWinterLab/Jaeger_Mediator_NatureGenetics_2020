########################################################################
# Publically available H3K27ac and H3K4me3 ChIP-seq data (GSE72622)
########################################################################

cd ${HCT116_working_dir}/chip_public

###########################################################################
# 1. Download raw data (H3K27ac, H3K4me3, input ChIP-seq on HCT116 cells)
###########################################################################

module load sratoolkit/2.8.2-1

# H3K27ac 
fastq-dump --split-files --gzip SRR2229187 SRR2229186 SRR2229185 SRR2229184

# H3K4me3
fastq-dump --outdir fastq_H3K4me3 --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR2229176 SRR2229177 SRR2229178 SRR2229179 

# Input
fastq-dump --outdir fastq_input --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR2229188 SRR2229189 SRR2229190 SRR2229191



###########################################################################
# 2. Processing the raw data (Mapping and filtering)
###########################################################################

module load bowtie/2.3.4
module load sambamba/0.5.5
module load htslib/1.8
module load samtools/1.8
module load bedtools/2.20.1 


genome=hg38



chip_data_name=HCT116_H3K27ac
for input_fastq in SRR2229184_pass.fastq.gz SRR2229185_pass.fastq.gz SRR2229186_pass.fastq.gz SRR2229187_pass.fastq.gz

#chip_data_name=HCT116_H3K4me3
#for input_fastq in SRR2229176_pass.fastq.gz SRR2229177_pass.fastq.gz SRR2229178_pass.fastq.gz SRR2229179_pass.fastq.gz

#chip_data_name=HCT116_input
#for input_fastq in SRR2229188_pass.fastq.gz SRR2229189_pass.fastq.gz SRR2229190_pass.fastq.gz SRR2229191_pass.fastq.gz

do

sample_name=$(echo ${chip_data_name}_${input_fastq%_pass.fastq.gz})

bowtie2 --very-sensitive --no-discordant -p 8 -x /data/groups/lab_bock/shared/resources/genomes/${genome}/indexed_bowtie2/${genome} --met-file ${sample_name}.aln_metrics.${genome}.txt -r ${input_fastq} 2> ${sample_name}.aln_rates.${genome}.txt | samtools view -S -b - | samtools sort -o ${sample_name}.trimmed.bowtie2.${genome}.bam -

sambamba markdup -t 8 -r --compression-level=0 ${sample_name}.trimmed.bowtie2.${genome}.bam ${sample_name}.trimmed.bowtie2.${genome}.filtered.nodups.nofilter.bam 2> ${sample_name}.dups_metrics.${genome}.txt
sambamba view -t 8 -f bam --valid -F "not unmapped and not (secondary_alignment or supplementary) and mapping_quality >= 30" ${sample_name}.trimmed.bowtie2.${genome}.filtered.nodups.nofilter.bam |sambamba sort -t 8 /dev/stdin -o ${sample_name}.trimmed.bowtie2.${genome}.filtered.bam

samtools index ${sample_name}.trimmed.bowtie2.${genome}.bam
samtools index ${sample_name}.trimmed.bowtie2.${genome}.filtered.bam

sambamba flagstat ${sample_name}.trimmed.bowtie2.${genome}.bam > ${sample_name}.trimmed.bowtie2.${genome}.flagstat
sambamba flagstat ${sample_name}.trimmed.bowtie2.${genome}.filtered.bam > ${sample_name}.trimmed.bowtie2.${genome}.filtered.flagstat


samtools sort ${sample_name}.trimmed.bowtie2.${genome}.filtered.bam > ${sample_name}.trimmed.bowtie2.${genome}.filtered.sorted.bam

done



# merge separate bam files
samtools merge ${chip_data_name}.filtered.bam *.filtered.sorted.bam
mv *.filtered.sorted.bam separate_files_${chip_data_name}

# index
samtools index ${chip_data_name}.filtered.bam

# generate bigWig
bedtools bamtobed -i ${chip_data_name}.filtered.bam | bedtools slop -i stdin -g /data/groups/lab_bock/shared/resources/genomes/${genome}/${genome}.chromSizes -s -l 0 -r 130 | /home/himrichova/1_Projects/0_scripts/fix_bedfile_genome_boundaries.py ${genome} | genomeCoverageBed -bg -g /data/groups/lab_bock/shared/resources/genomes/${genome}/${genome}.chromSizes -i stdin > ${chip_data_name}_${genome}.cov

awk 'NR==FNR{sum+= $4; next}{ $4 = ($4 / sum) * 1000000; print}' ${chip_data_name}_${genome}.cov ${chip_data_name}_${genome}.cov| sort -k1,1 -k2,2n > ${chip_data_name}_${genome}.normalized.cov

bedGraphToBigWig ${chip_data_name}_${genome}.normalized.cov /data/groups/lab_bock/shared/resources/genomes/${genome}/${genome}.chromSizes ${chip_data_name}_${genome}.bigWig

chmod 755 ${sample_name}_${genome}.bigWig






################################# 
# call peaks (hg38)
################################# 


module load MACS/2.1.0
#macs2 --version
#macs2 2.1.1.20160309

macs2 callpeak -t HCT116_H3K27ac.filtered.bam  -c HCT116_input.filtered.bam --fix-bimodal --extsize 180 --bw 200 --p 1e-09 -g hs -n HCT116_H3K27ac_hg38 --outdir ${HCT116_working_dir}/chip_public/bed

macs2 callpeak -t HCT116_H3K4me3.filtered.bam -c HCT116_input.filtered.bam --fix-bimodal --extsize 180 --bw 200 --p 1e-09 -g hs -n HCT116_H3K4me3_hg38 --outdir ${HCT116_working_dir}/chip_public/bed


awk -F "\t" '{{print $1 "\t" $4 "\t" "constituent_enhancer" "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $6 "\t" $4}}' HCT116_H3K27ac_hg38_peaks.narrowPeak > HCT116_H3K27ac_hg38_peaks.narrowPeak.log10_p9.gff
awk -F "\t" '{{print $1 "\t" $4 "\t" "constituent_enhancer" "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $6 "\t" $4}}' HCT116_H3K4me3_hg38_peaks.narrowPeak > HCT116_H3K4me3_hg38_peaks.narrowPeak.log10_p9.gff

