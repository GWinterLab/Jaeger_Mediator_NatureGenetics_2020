##########################################
# TT-seq data generated on KBM7 cells
##########################################

# Input: merged DMSO replicates - both 1h and 2h, only read1 (sense)

##############################################################################
# Map TTseq reads to genes (hg38) 
##############################################################################
# Coverage (bamToBed -> extended 100 bp -> coverage)


OUTPUT_DIR=coverage_TTseq_on_hg38_UCSC_full_length
REGIONS=hg38_UCSC_full_length.bed
#bedtools v2.20.1

mkdir ${OUTPUT_DIR}


for B in TTseq_KBM7_DMSO1h2h_read1_merged.bam
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



############################################################################
# Caluculate max RPKM per gene 
############################################################################

# scaling_factor = 1872.624 (total number of reads mapped to genes divided by 1,000,000)
awk -F "\t" '{print $0 "\t" $7/1872.624 "\t" $9/1000}' coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.bed |awk -F "\t" '{print $0 "\t" $11/$12}' > coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.bed



# Consider only protein coding genes (NM_)
multigrep.sh -f 5 -p "NM_" coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.bed > coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.bed


for gene_ID in $(cat coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.bed |cut -f4 |sort -u)
do
multigrep.sh -f 4 -w -p $gene_ID coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.bed |sort -k 13 -g -r |head -1 >> max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.bed
done


cut -f1-5,13 max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.bed > TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.max_RPKM_per_gene.only_NM.bed

#median expression value = 0.11015
awk -F "\t" '{if($13>=0.11015){print $0}}' max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.bed | sort -k13 -gr > max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.expressed_medianBased.bed 



cut -f4,13 max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.expressed_medianBased.bed|sort -k2 -gr > max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.expressed_medianBased.txt
awk -F "\t" '{print $1 "\t" $2 "\t" $2}' max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.expressed_medianBased.txt > max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.expressed_medianBased.3c.txt
