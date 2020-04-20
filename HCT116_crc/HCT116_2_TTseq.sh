##########################################
# TT-seq data generated on HCT116 cells
##########################################

# Input: merged DMSO replicates - both 1h and 2h, only read1 (sense)

samtools merge TTseq_HCT116_DMSO1h2h_read1_merged.bam L_with_x_nnnnnn_y_HCT116_DMSO1h_read1.bam L_with_x_nnnnnn_y_HCT116_DMSO2h_read1.bam

# Sort and index
samtools sort TTseq_HCT116_DMSO1h2h_read1_merged.bam > TTseq_HCT116_DMSO1h2h_read1_merged.sorted.bam

samtools index TTseq_HCT116_DMSO1h2h_read1_merged.sorted.bam





##############################################################################
# Map TTseq reads to genes (hg38)
##############################################################################
# Coverage (bamToBed -> extended 100 bp -> coverage)

OUTPUT_DIR=coverage_HCT116_TTseq_1h2h_on_hg38_UCSC_full_length
REGIONS=hg38_UCSC_full_length.bed
#bedtools v2.20.1

mkdir ${OUTPUT_DIR}


for B in TTseq_HCT116_DMSO1h2h_read1_merged.sorted.bam
do
echo "${B}"

bedtools bamtobed -i ${B} | bedtools slop -g /data/groups/lab_bock/shared/resources/genomes/hg38/hg38.chromSizes -b 100 > ${OUTPUT_DIR}/$(basename ${B%.bam}).extended100bp.bed

if [ ! -f ${OUTPUT_DIR}/coverage-$(basename ${B%.bam}).extended100bp.bed ]
then
echo "coverage-$(basename ${B%.bam}).extended100bp.bed does NOT exist yet."
coverageBed -a ${OUTPUT_DIR}/$(basename ${B%.bam}).extended100bp.bed -b ${REGIONS} > ${OUTPUT_DIR}/coverage-$(basename ${B%.bam}).extended100bp.bed
else
echo "$(basename ${B%.bam}).extended100bp.bed exists."
fi
done


######################################
# Caluculate RPKM
######################################

# scaling_factor = 865.682 (total number of reads mapped to genes divided by 1,000,000)


cat ${coverage_file}|awk -F "\t" '{print $0 "\t" $7/865.682 "\t" $8/1000}' |awk -F "\t" '{print $0 "\t" $9/$10}' > ${coverage_file%.bed}.RPKM.bed




######################################
# max RPKM per gene:
######################################


# only protein coding NM_
coverage_RPKM_file=cov_TTseq_HCT116_DMSO1h2h_read1_merged.extended100bp.RPKM.bed


multigrep.sh -f 5 -p "NM_"  ${coverage_RPKM_file} > ${coverage_RPKM_file}.only_NM.bed


for gene_ID in $(cat ${coverage_RPKM_file}.only_NM.bed |cut -f4 |sort -u)
do
#echo $gene_ID
multigrep.sh -f 4 -w -p $gene_ID ${coverage_RPKM_file}.only_NM.bed |sort -k 11 -g -r |head -1 >> max_RPKM_per_gene.${coverage_RPKM_file}.only_NM.bed
done






cut -f1-5,11 max_RPKM_per_gene.${coverage_RPKM_file}.only_NM.bed > TTseq_${coverage_RPKM_file}.max_RPKM_per_gene.only_NM.bed


# identify median
# from the max.expressed transcripts per gene
# median expression value = 0.108945

cat max_RPKM_per_gene.${coverage_RPKM_file}.only_NM.bed | awk -F "\t" '{if($11>=0.108945){print $0}}' | sort -k11 -gr > max_RPKM_per_gene.${coverage_RPKM_file}.only_NM.expressed_medianBased.bed 



less max_RPKM_per_gene.${coverage_RPKM_file}.only_NM.expressed_medianBased.bed|cut -f4,11|sort -k2 -gr > max_RPKM_per_gene.${coverage_RPKM_file}.only_NM.expressed_medianBased.txt

awk -F "\t" '{print $1 "\t" $2 "\t" $2}' max_RPKM_per_gene.${coverage_RPKM_file}.only_NM.expressed_medianBased.txt > max_RPKM_per_gene.${coverage_RPKM_file}.only_NM.expressed_medianBased.3c.txt
