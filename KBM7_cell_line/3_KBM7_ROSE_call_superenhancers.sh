#####################################################
# Generate super-enhancers using the Rose algorithm
#####################################################



#######################################################################################################
# 1. Remove the constitutive enhancers that overlap with H3K4me3&TSS (i.e. remove potential promoters)
#######################################################################################################

#intersectBed
#Version: v2.20.1

TSS_bed=hg38_UCSC.TSS.sorted.bed

intersectBed -a KBM7_H3K27ac_hg38_peaks.narrowPeak.log10_p9.bed -b <(intersectBed -a KBM7_H3K4me3_hg38_peaks.narrowPeak.log10_p9.bed -b <(bedtools slop -i ${TSS_bed} -g /data/groups/lab_bock/shared/resources/genomes/hg38/hg38.chromSizes -b 1000 ) -wa |sort -u) -wa -v |sort -u > KBM7_H3K27ac_hg38_peaks.noOverlap_H3K4me3_TSS.bed 

awk -F "\t" '{if($8>=9){print $1 "\t" $4 "\t" "constituent_enhancer" "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $6 "\t" $4}}' KBM7_H3K27ac_hg38_peaks.noOverlap_H3K4me3_TSS.bed > KBM7_H3K27ac_hg38_peaks.noOverlap_H3K4me3_TSS.gff


##########################################################
# 2. Call super-enhancers using the Rose tool (hg38)
##########################################################
# download HG38 (from https://raw.githubusercontent.com/linlabbcm/rose2/master/rose2/annotation/hg38_refseq.ucsc) into dir ~/software/young_computation-rose-1a9bb86b5464_adapted/annotation 

cd ~/software/young_computation-rose-1a9bb86b5464_adapted


enhancers=${KBM7_working_dir}/bed/KBM7_H3K27ac_hg38_peaks.noOverlap_H3K4me3_TSS.gff
K27ac_bam=${KBM7_working_dir}/mapping/KBM7_H3K27ac.trimmed.bowtie2.hg38.filtered.bam
input_bam=${KBM7_working_dir}/mapping/KBM7_input.trimmed.bowtie2.hg38.filtered.bam
outfolder=${KBM7_working_dir}/bed/Rose_SE_hg38_p9_noOverlap_H3K4me3_TSS/

python ROSE_main.py -c $input_bam -g HG38 -i $enhancers -r $K27ac_bam -o $outfolder


