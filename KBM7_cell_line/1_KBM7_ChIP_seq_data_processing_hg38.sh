################################# 
# Map chip-seq data to hg38
################################# 

${chip_scripts}/1_chip_seq_mapping.sh hg38 KBM7_input.trimmed.fastq KBM7_input

${chip_scripts}/1_chip_seq_mapping.sh hg38 KBM7_H3K4me3.trimmed.fastq KBM7_H3K4me3

${chip_scripts}/1_chip_seq_mapping.sh hg38 KBM7_H3K27ac.trimmed.fastq KBM7_H3K27ac


################################# 
# Call peaks (hg38)
################################# 


module load MACS/2.1.0
#macs2 --version
#macs2 2.1.1.20160309

macs2 callpeak -t KBM7_Hq3K27ac.trimmed.bowtie2.hg38.filtered.bam -c KBM7_input.trimmed.bowtie2.hg38.filtered.bam --fix-bimodal --extsize 180 --bw 200 --p 1e-09 -g hs -n KBM7_H3K27ac_hg38 --outdir bed

macs2 callpeak -t KBM7_H3K4me3.trimmed.bowtie2.hg38.filtered.bam -c KBM7_input.trimmed.bowtie2.hg38.filtered.bam --fix-bimodal --extsize 180 --bw 200 --p 1e-09 -g hs -n KBM7_H3K4me3_hg38 --outdir bed


# create gff files:
awk -F "\t" '{{print $1 "\t" $4 "\t" "constituent_enhancer" "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $6 "\t" $4}}' KBM7_H3K27ac_hg38_peaks.narrowPeak.log10_p9.bed > KBM7_H3K27ac_hg38_peaks.narrowPeak.log10_p9.gff
awk -F "\t" '{{print $1 "\t" $4 "\t" "constituent_enhancer" "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $6 "\t" $4}}' KBM7_H3K27ac_hg38_peaks.narrowPeak.log10_p9.bed > KBM7_H3K4me3_hg38_peaks.narrowPeak.log10_p9.gff

