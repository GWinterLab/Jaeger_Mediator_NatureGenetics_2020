########################################################
# Identification of autoregulatory TFs for HCT116
########################################################

############################
# TADs HCT116

TADs_bed=HCT-116_RAD21-mAC_no_auxin.Rao_2017-raw.domains.adapted.bed

#####################################################################################################################################################################
# Expressed genes defined using TT-seq HCT116 - merged DMSO_1h and DMSO_2h (only read1), median based - from the max.expression per gene (protein coding genes)

expressed_genes_bed=expressed_TTseq_onlyNM.bed

####################################################################################
# Super-enhancers and subpbeaks defined using public Mol Cell ChIP-seq data

SE_bed=HCT116_H3K27ac_hg38_peaks_Gateway_SuperEnhancers.bed

subpeaks=HCT116_H3K27ac_hg38_peaks.noOverlap_H3K4me3_TSS.bed



awk -F "\t" '{print $0 "\t" $1":"$2"-"$3}' ${SE_bed} > SEs_faIDs.bed




##############################################################################
# 1. Assign expressed genes to SEs
##############################################################################

#bedtools v2.20.1

# generate a file including SEs that are in the same TAD as an expressed gene:
bedtools intersect -a ${TADs_bed} -b ${SE_bed} -wo > TADs_with_SE.bed

bedtools intersect -a TADs_with_SE.bed -b ${expressed_genes_bed} -wo|cut -f4-|sort -u > TADs_with_SE_and_expressed_genes.bed



##############################################################################
# 2. Make fasta file of extended enhancers in SEs (hg38)
##############################################################################


intersectBed -a <(bedtools slop -i ${subpeaks} -g /data/groups/lab_bock/shared/resources/genomes/hg38/hg38.chromSizes -b 500) -b TADs_with_SE_and_expressed_genes.bed -wa |sort -u > subpeaks_500bp_extended_in_SEs_assigned_to_expressed_genes.bed


bedtools getfasta -fi /data/groups/lab_bock/shared/resources/genomes/hg38/hg38.fa -bed subpeaks_500bp_extended_in_SEs_assigned_to_expressed_genes.bed -fo subpeaks_500bp_extended_in_SEs_assigned_to_expressed_genes.fasta








##############################################################################
# 3. Select motifs of expressed TFs
##############################################################################

# Motif database used by iRegulon (Janky et al., 2014) and i-cisTarget (Imrichova et al., 2015), where TFs identified as "Unlikely to be sequence specific TF" or "ssDNA/RNA binding" are filetered out (Lambert et al., 2018).
# motifs-v8-nr.hgnc-m0.000-o0.0.directly.true_TFs.tbl

cat motifs-v8-nr.hgnc-m0.000-o0.0.directly.true_TFs.tbl |multigrep.sh -f 6 -w -g <(cut -f4 ${expressed_genes_bed}) |cut -f6 |sort -u > expressed_TFs_list.txt


# Assign expressed TFs to SEs:
multigrep.sh -f 11 -w -g expressed_TFs_list.txt TADs_with_SE_and_expressed_genes.bed |cut -f11 |sort -u > expressed_TFs_list.assigned_to_SEs.txt

multigrep.sh -f 6 -w -g expressed_TFs_list.assigned_to_SEs.txt motifs-v8-nr.hgnc-m0.000-o0.0.directly.true_TFs.tbl |cut -f1 |sort -u > motifs-v8-nr.hgnc-m0.000-o0.0.directly.expressed.assigned.txt







#####################################################################################################################################
# 4. score all SEs with the motifs annotated for expressed TFs that are in the same TAD with any SE (cbust)
#####################################################################################################################################
# Version: CLUSTER-BUSTER:  2018-07-04 19:09:53 -0400  commit: 0ee2a95

#module load cluster-buster/2018

mkdir cbust_output
for motif_name in $(cat motifs-v8-nr.hgnc-m0.000-o0.0.directly.expressed.assigned.txt)
do
cb_file=${Motifs}/motifCollection_v8/motif_collection_v8/singletons/${motif_name}.cb

echo $motif_name
cbust -c 5 -m 6 -G 0 -f 5 ${cb_file} subpeaks_500bp_extended_in_SEs_assigned_to_expressed_genes.fasta > cbust_output/SEs.${motif_name}.c5m6G0f5.cbust2018_out.bed
done



###################################################################
# 5. extract cbust results on extended Ac peaks
###################################################################
mkdir cbust_output_SE_IDs

SE_IDs=SEs_faIDs.bed

for i in cbust_output/*bed
do
echo $(basename ${i})

bedtools intersect -a cbust_output/$(basename ${i}) -b ${SE_IDs} -wo > cbust_output_SE_IDs/$(basename ${i%.bed}).on_extended_peaks.bed
done




###################################################################################################
# 6. HCT116 AUTO TFs
#
# Identify Auto-regulatory TFs (TFs that have their own TF-binding motif in their assigned SE)
###################################################################################################

cbust_on_SEs_directory=${HCT116_working_dir}/chip_public/bed/Rose_SE_hg38_p9_noOverlap_H3K4me3_TSS/CRC_inhouse/cbust_output_SE_IDs/
TFs_list=${HCT116_working_dir}/chip_public/bed/Rose_SE_hg38_p9_noOverlap_H3K4me3_TSS/CRC_inhouse/expressed_TFs_list.assigned_to_SEs.txt


directly_annot_motifs_tbl=${Motifs}/cleaned_motif_database/motifs-v8-nr.hgnc-m0.000-o0.0.directly.true_TFs.tbl


awk -F "\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $1 ":" $2 "-" $3 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' TADs_with_SE_and_expressed_genes.bed > TADs_with_SE_and_TSS.bed



for TF in $(cat TADs_with_SE_and_TSS.bed | multigrep.sh -f 12 -w -g ${TFs_list}|cut -f12|sort -u)
do
rm TF_motifs.tmp
rm TF_SEs.tmp

cat TADs_with_SE_and_TSS.bed | multigrep.sh -f 12 -w -p ${TF} |cut -f1-7 |sort -u > TF_SEs.tmp
cat ${directly_annot_motifs_tbl} | multigrep.sh -f 6 -w -p ${TF} |cut -f1 |sort -u > TF_motifs.tmp

for SE_name in $(cat TF_SEs.tmp |cut -f7 |sort -u)
do
echo $SE_name

for motif_ID in $(cat TF_motifs.tmp)
do
echo $motif_ID

paste <(echo $TF) <(echo $SE_name) <(echo $motif_ID) <(multigrep.sh -f 23 -w -p ${SE_name} ${cbust_on_SEs_directory}SEs.${motif_ID}.c5m6G0f5.cbust2018_out.on_extended_peaks.bed | multigrep.sh -f 11 -w -p "motif"|cut -f1-4|sort -u|wc -l) <(multigrep.sh -f 23 -w -p ${SE_name} ${cbust_on_SEs_directory}SEs.${motif_ID}.c5m6G0f5.cbust2018_out.on_extended_peaks.bed | multigrep.sh -f 11 -w -p "cluster"|cut -f1-3|sort -u|wc -l) >> auto_TF_results.TF_SE_motifID_motifs_clusters.txt

done
done
done





###########################
# get auto-regulated TFs based on specific threshold (number of motifs, number of motif clusters in the assigned SE)

awk -F "\t" '{if($4>=3 && $5>=3){print $0}}' auto_TF_results.TF_SE_motifID_motifs_clusters.txt > candidate_autoTFs.3motifs_3clusters.tbl

awk -F "\t" '{if($4>=6 && $5>=6){print $0}}' auto_TF_results.TF_SE_motifID_motifs_clusters.txt > candidate_autoTFs.6motifs_6clusters.tbl



##########################################################################
# 7. CRCs
# Core transcriptional regulatory circuitries (CRCs) inferring
# only closest genes
##########################################################################

cbust_on_SEs_directory=${HCT116_working_dir}/chip_public/bed/Rose_SE_hg38_p9_noOverlap_H3K4me3_TSS/CRC_inhouse/cbust_output_SE_IDs/

file_with_closest_genes=${HCT116_working_dir}/chip_public/bed/Rose_SE_hg38_p9_noOverlap_H3K4me3_TSS/CRC_inhouse/genesets/geneset_2b.txt


number_of_motifs=3
number_of_clusters=3

mkdir autoTFs_${number_of_motifs}motifs_${number_of_clusters}clusters_closest_genes_test
cd autoTFs_${number_of_motifs}motifs_${number_of_clusters}clusters_closest_genes_test

awk -F "\t" -v number_of_motifs=${number_of_motifs} -v number_of_clusters=${number_of_clusters} '{if($4>=number_of_motifs && $5>=number_of_clusters){print $0}}' ../auto_TF_results.TF_SE_motifID_motifs_clusters.txt |multigrep.sh -f 1 -w -g ${file_with_closest_genes} > candidate_autoTFs.${number_of_motifs}motifs_${number_of_clusters}clusters.tbl

table_of_candidate_autoregulatedTFs=candidate_autoTFs.${number_of_motifs}motifs_${number_of_clusters}clusters.tbl

rm candidate_SEs.txtCRC_results_motif_occurrence.txt
cut -f2 ${table_of_candidate_autoregulatedTFs} |sort -u > candidate_SEs.txt



for TF in $(cat ${table_of_candidate_autoregulatedTFs} |cut -f1 |sort -u)
do
rm TF_motifs.tmp
cat ${table_of_candidate_autoregulatedTFs} | multigrep.sh -f 1 -w -p ${TF} |cut -f3 |sort -u > TF_motifs.tmp

for SE_name in $(cat candidate_SEs.txt)
do

for motif_ID in $(cat TF_motifs.tmp)
do

for TG in $( cat $table_of_candidate_autoregulatedTFs | multigrep.sh -f 2 -w -p ${SE_name} |cut -f1|sort -u)
do 

paste <(echo ${TF}) <(echo $motif_ID) <(echo $SE_name) <(echo $TG) <(multigrep.sh -f 23 -w -p ${SE_name} ${cbust_on_SEs_directory}SEs.${motif_ID}.c5m6G0f5.cbust2018_out.on_extended_peaks.bed | multigrep.sh -f 11 -w -p "motif"|cut -f1-4|sort -u|wc -l) <(multigrep.sh -f 23 -w -p ${SE_name} ${cbust_on_SEs_directory}SEs.${motif_ID}.c5m6G0f5.cbust2018_out.on_extended_peaks.bed | multigrep.sh -f 11 -w -p "cluster"|cut -f1-3|sort -u|wc -l) >> CRC_results_motif_occurrence.txt

done
done
done
done


cp CRC_results_motif_occurrence.txt CRC_results_motif_occurrence.motif${number_of_motifs}_cluster${number_of_clusters}.txt



awk -F "\t" -v number_of_motifs=${number_of_motifs} -v number_of_clusters=${number_of_clusters} '{if($5>=number_of_motifs && $6>=number_of_clusters){print $0}}' CRC_results_motif_occurrence.txt > CRC_results_motif_occurrence.m_c.txt




