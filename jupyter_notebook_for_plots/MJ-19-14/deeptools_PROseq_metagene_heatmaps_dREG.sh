################################################################
##### dREG was run on MJ-19-14_DMSO_1h_merged bigWig files #####
################################################################

############################################
##### get normal and SE gene BED files #####
############################################

NM_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.bed
#NR_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.bed
SE_list=/scratch/lab_winter/martin/Mediator_CRC/SEs_464_TTseq_NMonly.txt
autoreg_list=/scratch/lab_winter/martin/Mediator_CRC/CRC_24_TTseq_NMonly_autoregulatory_6motifs.txt

SE_BED=/scratch/lab_winter/martin/Mediator_CRC/BED_for_deeptools/SEs_464_TTseq_NMonly.bed
other_BED=/scratch/lab_winter/martin/Mediator_CRC/BED_for_deeptools/not_SEs_9276_TTseq_NMonly.bed
autoreg_BED=/scratch/lab_winter/martin/Mediator_CRC/BED_for_deeptools/CRC_24_TTseq_NMonly_autoregulatory_6motifs.bed

awk '(NR==FNR) {genes[$0]; next}; ($4 in genes) && ($13 > 0.110155) {print $0}' $SE_list $NM_cov | cut -f1-6 >  $SE_BED #cutoff = median TT-seq RPKM
awk '(NR==FNR) {genes[$0]; next}; !($4 in genes) && ($13 > 0.110155) {print $0}' $SE_list $NM_cov | cut -f1-6 >  $other_BED
awk '(NR==FNR) {genes[$0]; next}; ($4 in genes) && ($13 > 0.110155) {print $0}' $autoreg_list $NM_cov | cut -f1-6 >  $autoreg_BED #cutoff = median TT-seq RPKM


#################################################################
##### use deeptools 3.3.0 to make the matrices and heatmaps #####
#################################################################
DMSO_1h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/bigWig/MJ-19-14_DMSO_1h_merged_plus_dm6_normalized.bw
DMSO_1h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/bigWig/MJ-19-14_DMSO_1h_merged_minus_dm6_normalized.bw
dTAG_1h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/bigWig/MJ-19-14_dTAG7_1h_merged_plus_dm6_normalized.bw
dTAG_1h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/bigWig/MJ-19-14_dTAG7_1h_merged_minus_dm6_normalized.bw

## computeMatrix for the plus and minus read signals separately
plus_reads_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_plus_reads_autoreg.matrix.gz
minus_reads_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_minus_reads_autoreg.matrix.gz

computeMatrix scale-regions -S $DMSO_1h_plus_bw $dTAG_1h_plus_bw -R $autoreg_BED $SE_BED $other_BED -o $plus_reads_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel DMSO_1h_plus dTAG7_1h_plus -p max/2 -bs 25
computeMatrix scale-regions -S $DMSO_1h_minus_bw $dTAG_1h_minus_bw -R $autoreg $SE_BED $other_BED -o $minus_reads_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel DMSO_1h_minus dTAG7_1h_minus -p max/2 --scale "-1" -bs 25

## retain only signal from the respective gene's strand (and also antisense)
plus_strand_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_plus_strand_autoreg.matrix.gz
minus_strand_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_minus_strand_autoreg.matrix.gz
plus_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_plus_antisense_autoreg.matrix.gz
minus_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_minus_antisense_autoreg.matrix.gz

computeMatrixOperations filterStrand -m $plus_reads_matrix -o $plus_strand_matrix --strand +
computeMatrixOperations filterStrand -m $minus_reads_matrix -o $minus_strand_matrix --strand -
computeMatrixOperations filterStrand -m $plus_reads_matrix -o $plus_antisense_matrix --strand -
computeMatrixOperations filterStrand -m $minus_reads_matrix -o $minus_antisense_matrix --strand +

## merge the strands back together to get the "sense signal" matrix
merged_sense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_sense_autoreg.matrix.gz
merged_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_antisense_autoreg.matrix.gz

computeMatrixOperations rbind -m $plus_strand_matrix $minus_strand_matrix -o $merged_sense_matrix
computeMatrixOperations rbind -m $plus_antisense_matrix $minus_antisense_matrix -o $merged_antisense_matrix

## finally plot the sense PRO-seq signal heatmaps & metagenes
sense_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_sense_heatmap_autoreg.pdf
antisense_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_antisense_heatmap_autoreg.pdf

plotHeatmap -m $merged_sense_matrix -o $sense_output --dpi 300 --colorList "#ffffff00,#960505" --alpha 1 --regionsLabel "24 autoreg TFs" "464 SEs" "9276 non-SEs" --sortUsing mean --sortUsingSamples 1 #--zMax 0.2 --zMin "-0.2"
plotHeatmap -m $merged_antisense_matrix -o $antisense_output --dpi 300 --colorList "#ffffff00,#050596" --alpha 1 --regionsLabel "24 autoreg TFs" "464 SEs" "9276 non-SEs" --sortUsing mean --sortUsingSamples 1 #--zMax 0.2 --zMin "-0.2"




##########################################################
##### use deeptools 3.3.0 to make dTAG/DMSO profiles #####
##########################################################
sense_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_sense_profile_fill_autoreg.pdf
antisense_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_antisense_profile_fill_autoreg.pdf

plotProfile -m $merged_sense_matrix -o $sense_profile_output --dpi 300 --colors "#faa41a" "#be1e2d" "#231f20" --plotType fill
plotProfile -m $merged_antisense_matrix -o $antisense_profile_output --dpi 300 --colors "#faa41a" "#be1e2d" "#231f20" --plotType fill

## now seperate by SE vs. non-SE rather than condition
sense_perGroup_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_sense_perGroup_profile_fill_autoreg.pdf
antisense_perGroup_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_antisense_perGroup_profile_fill_autoreg.pdf

plotProfile -m $merged_sense_matrix -o $sense_perGroup_profile_output --dpi 300 --colors "#be1e2d" "#231f20" --plotType fill --perGroup
plotProfile -m $merged_antisense_matrix -o $antisense_perGroup_profile_output --dpi 300 --colors "#be1e2d" "#231f20" --plotType fill --perGroup

# raw values for manual plotting
plotProfile -m $merged_sense_matrix -o $sense_perGroup_profile_output --dpi 300 --plotWidth 10 --plotHeight 10 --colors "#be1e2d" "#231f20" --plotType fill --perGroup --outFileNameData MJ-19-14_metagenes_autoregTF_6cell-smoothed_profile_values_raw_output.txt



#################################################################
##### use deeptools 3.3.0 to make log2FC dTAG/DMSO heatmaps #####
#################################################################
BAM_1h=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/aligned/MJ-19-14_dTAG7_1h_merged_1bp.sorted.bam
BAM_DMSO=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/aligned/MJ-19-14_DMSO_1h_merged_1bp.sorted.bam
scale_1h=0.176419
scale_DMSO=0.247928

log2FC_1h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/bigWig/MJ-19-14_log2FC_1h-over-DMSO_plus_normalized.bw
log2FC_1h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/bigWig/MJ-19-14_log2FC_1h-over-DMSO_minus_normalized.bw

## compute log2FC bigwig files separately per strand
bamCompare -b1 $BAM_1h -b2 $BAM_DMSO -o $log2FC_1h_plus_bw --scaleFactors $scale_1h:$scale_DMSO --operation log2 --pseudocount 0.01 -bs 1 -p max/2 --samFlagExclude 16
bamCompare -b1 $BAM_1h -b2 $BAM_DMSO -o $log2FC_1h_minus_bw --scaleFactors $scale_1h:$scale_DMSO --operation log2 --pseudocount 0.01 -bs 1 -p max/2 --samFlagInclude 16

## computeMatrix for the plus and minus read signals separately
plus_reads_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_plus_reads_log2FC.matrix.gz
minus_reads_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_minus_reads_log2FC.matrix.gz

computeMatrix scale-regions -S $log2FC_1h_plus_bw -R $SE_BED $other_BED -o $plus_reads_log2_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel log2FC_1h_vs_DMSO_plus -p max/2 -bs 25
computeMatrix scale-regions -S $log2FC_1h_minus_bw -R $SE_BED $other_BED -o $minus_reads_log2_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel log2FC_1h_vs_DMSO_minus -p max/2 -bs 25

## retain only signal from the respective gene's strand (and also antisense)
plus_strand_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_plus_strand_log2FC.matrix.gz
minus_strand_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_minus_strand_log2FC.matrix.gz
plus_antisense_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_plus_antisense_log2FC.matrix.gz
minus_antisense_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_minus_antisense_log2FC.matrix.gz

computeMatrixOperations filterStrand -m $plus_reads_log2_matrix -o $plus_strand_log2_matrix --strand +
computeMatrixOperations filterStrand -m $minus_reads_log2_matrix -o $minus_strand_log2_matrix --strand -
computeMatrixOperations filterStrand -m $plus_reads_log2_matrix -o $plus_antisense_log2_matrix --strand -
computeMatrixOperations filterStrand -m $minus_reads_log2_matrix -o $minus_antisense_log2_matrix --strand +

## merge the strands back together to get the "sense signal" matrix
merged_sense_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_sense_log2FC.matrix.gz
merged_antisense_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_antisense_log2FC.matrix.gz

computeMatrixOperations rbind -m $plus_strand_log2_matrix $minus_strand_log2_matrix -o $merged_sense_log2_matrix
computeMatrixOperations rbind -m $plus_antisense_log2_matrix $minus_antisense_log2_matrix -o $merged_antisense_log2_matrix

## finally plot the sense PRO-seq signal heatmaps & metagenes
sense_log2_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_sense_log2FC_heatmap.pdf
antisense_log2_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_antisense_log2FC_heatmap.pdf
both_log2_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_both_log2FC_heatmap.pdf

plotHeatmap -m $merged_sense_log2_matrix -o $sense_log2_output --dpi 300 --colorList "#050596,#ffffff00,#960505" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortRegions ascend --zMax 1 --zMin "-1"
plotHeatmap -m $merged_antisense_log2_matrix -o $antisense_log2_output --dpi 300 --colorList "#050596,#ffffff00,#960505" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortRegions ascend --zMax 1 --zMin "-1"





################################################
##### combined DMSO | 1h | log2FC heatmaps #####
################################################
## merge left-to-right sample signal and log2FC
sense_master_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_sense_DMSO-1h-log2FC.matrix.gz

computeMatrixOperations cbind -m $merged_sense_matrix $merged_sense_log2_matrix -o $sense_master_matrix

## plot combined DMSO | 1h | log2FC signal heatmaps
sense_master_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_metagene_merged_sense_DMSO-1h-log2FC_heatmap.pdf

plotHeatmap -m $sense_master_matrix -o $sense_master_output --dpi 300 --colorList "#ffffff00,#960505" "#ffffff00,#960505" "#f2de05,#231f20" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortUsingSamples 1 --zMax 0.2 0.2 0 --zMin 0 0 "-1"





#################################################################
##### similar heatmaps/metagenes for NVP2 combo experiments #####
#################################################################
#################################################################
##### use deeptools 3.3.0 to make the matrices and heatmaps #####
#################################################################

DMSO_2h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_DMSO_2h_merged_plus_dm6_normalized.bw
DMSO_2h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_DMSO_2h_merged_minus_dm6_normalized.bw
NVP2_30min_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/bigWig/MJ-19-30_NVP2_30min_merged_plus_dm6_normalized.bw
NVP2_30min_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/bigWig/MJ-19-30_NVP2_30min_merged_minus_dm6_normalized.bw
dTAG_2h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/bigWig/MJ-19-30_dTAG7_2h_merged_plus_dm6_normalized.bw
dTAG_2h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/bigWig/MJ-19-30_dTAG7_2h_merged_minus_dm6_normalized.bw
dTAG_NVP2_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/bigWig/MJ-19-30_dTAG7_NPV2_combo_merged_plus_dm6_normalized.bw
dTAG_NVP2_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/bigWig/MJ-19-30_dTAG7_NPV2_combo_merged_minus_dm6_normalized.bw

## computeMatrix for the plus and minus read signals separately
plus_reads_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_plus_reads.matrix.gz
minus_reads_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_minus_reads.matrix.gz

computeMatrix scale-regions -S $DMSO_2h_plus_bw $NVP2_30min_plus_bw $dTAG_2h_plus_bw $dTAG_NVP2_plus_bw -R $SE_BED $other_BED -o $plus_reads_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel MJ-19-30_DMSO_2h_plus MJ-19-30_NVP2_30min_plus MJ-19-30_dTAG7_2h_plus MJ-19-30_dTAG7_NVP2_combined_plus -p max/2 -bs 25
computeMatrix scale-regions -S $DMSO_2h_minus_bw $NVP2_30min_minus_bw $dTAG_2h_minus_bw $dTAG_NVP2_minus_bw -R $SE_BED $other_BED -o $minus_reads_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel MJ-19-30_DMSO_2h_minus MJ-19-30_NVP2_30min_minus MJ-19-30_dTAG7_2h_minus MJ-19-30_dTAG7_NVP2_combined_minus -p max/2 --scale "-1" -bs 25

## retain only signal from the respective gene's strand (and also antisense)
plus_strand_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_plus_strand.matrix.gz
minus_strand_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_minus_strand.matrix.gz
plus_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_plus_antisense.matrix.gz
minus_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_minus_antisense.matrix.gz

computeMatrixOperations filterStrand -m $plus_reads_matrix -o $plus_strand_matrix --strand +
computeMatrixOperations filterStrand -m $minus_reads_matrix -o $minus_strand_matrix --strand -
computeMatrixOperations filterStrand -m $plus_reads_matrix -o $plus_antisense_matrix --strand -
computeMatrixOperations filterStrand -m $minus_reads_matrix -o $minus_antisense_matrix --strand +

## merge the strands back together to get the "sense signal" matrix
merged_sense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_merged_sense.matrix.gz
merged_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_merged_antisense.matrix.gz

computeMatrixOperations rbind -m $plus_strand_matrix $minus_strand_matrix -o $merged_sense_matrix
computeMatrixOperations rbind -m $plus_antisense_matrix $minus_antisense_matrix -o $merged_antisense_matrix

## finally plot the sense PRO-seq signal heatmaps & metagenes
sense_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_merged_sense_heatmap.pdf
antisense_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_merged_antisense_heatmap.pdf

plotHeatmap -m $merged_sense_matrix -o $sense_output --dpi 300 --colorList "#ffffff00,#960505" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortUsingSamples 1 --outFileSortedRegions sense_regions_4samples.bed #--zMax 0.2 --zMin "-0.2"
plotHeatmap -m $merged_antisense_matrix -o $antisense_output --dpi 300 --colorList "#ffffff00,#050596" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortUsingSamples 1 #--zMax 0.2 --zMin "-0.2"




#######################################
##### MJ-19-30 dTAG/DMSO profiles #####
#######################################
sense_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_merged_sense_profile_fill.pdf
antisense_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_merged_antisense_profile_fill.pdf

plotProfile -m $merged_sense_matrix -o $sense_profile_output --dpi 300 --colors "#be1e2d" "#231f20" --plotType fill
plotProfile -m $merged_antisense_matrix -o $antisense_profile_output --dpi 300 --colors "#be1e2d" "#231f20" --plotType fill

## now seperate by SE vs. non-SE rather than condition
sense_perGroup_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_merged_sense_perGroup_profile_fill.pdf
antisense_perGroup_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_merged_antisense_perGroup_profile_fill.pdf

plotProfile -m $merged_sense_matrix -o $sense_perGroup_profile_output --dpi 300 --colors "#231f2080" "#94919180" "#ed670780" "#be1e2d80" --plotType fill --perGroup
plotProfile -m $merged_antisense_matrix -o $antisense_perGroup_profile_output --dpi 300 --colors "#231f2080" "#94919180" "#ed670780" "#be1e2d80" --plotType fill --perGroup





#################################################################
##### use deeptools 3.3.0 to make log2FC dTAG/DMSO heatmaps #####
#################################################################
BAM_DMSO=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/aligned/MJ-19-30_DMSO_2h_merged_1bp.sorted.bam
BAM_NVP2=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/aligned/MJ-19-30_NVP2_30min_merged_1bp.sorted.bam
BAM_dTAG=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/aligned/MJ-19-30_dTAG7_2h_merged_1bp.sorted.bam
BAM_combo=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/aligned/MJ-19-30_dTAG7_NPV2_combo_merged_1bp.sorted.bam
scale_DMSO=0.412471
scale_NVP2=0.392713
scale_dTAG=0.274846
scale_combo=0.37138

log2FC_NVP2_vs_DMSO_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_log2FC_NVP2-over-DMSO_plus_normalized.bw
log2FC_NVP2_vs_DMSO_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_log2FC_NVP2-over-DMSO_minus_normalized.bw
log2FC_combo_vs_dTAG_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_log2FC_combo_vs_dTAG_plus_normalized.bw
log2FC_combo_vs_dTAG_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_log2FC_combo_vs_dTAG_minus_normalized.bw

## compute log2FC bigwig files separately per strand
bamCompare -b1 $BAM_NVP2 -b2 $BAM_DMSO -o $log2FC_NVP2_vs_DMSO_plus_bw --scaleFactors $scale_NVP2:$scale_DMSO --operation log2 --pseudocount 0.01 -bs 1 -p max/2 --samFlagExclude 16
bamCompare -b1 $BAM_NVP2 -b2 $BAM_DMSO -o $log2FC_NVP2_vs_DMSO_minus_bw --scaleFactors $scale_NVP2:$scale_DMSO --operation log2 --pseudocount 0.01 -bs 1 -p max/2 --samFlagInclude 16
bamCompare -b1 $BAM_combo -b2 $BAM_dTAG -o $log2FC_combo_vs_dTAG_plus_bw --scaleFactors $scale_combo:$scale_dTAG --operation log2 --pseudocount 0.01 -bs 1 -p max/2 --samFlagExclude 16
bamCompare -b1 $BAM_combo -b2 $BAM_dTAG -o $log2FC_combo_vs_dTAG_minus_bw --scaleFactors $scale_combo:$scale_dTAG --operation log2 --pseudocount 0.01 -bs 1 -p max/2 --samFlagInclude 16

## computeMatrix for the plus and minus read signals separately
plus_reads_NVP2_vs_DMSO_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagenes_plus_reads_log2FC.matrix.gz
minus_reads_NVP2_vs_DMSO_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagene_minus_reads_log2FC.matrix.gz
plus_reads_combo_vs_dTAG_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_plus_reads_log2FC.matrix.gz
minus_reads_combo_vs_dTAG_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagene_minus_reads_log2FC.matrix.gz

computeMatrix scale-regions -S $log2FC_NVP2_vs_DMSO_plus_bw -R $SE_BED $other_BED -o $plus_reads_NVP2_vs_DMSO_log2_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel log2FC_NVP2_vs_DMSO_plus -p max/2 -bs 25
computeMatrix scale-regions -S $log2FC_NVP2_vs_DMSO_minus_bw -R $SE_BED $other_BED -o $minus_reads_NVP2_vs_DMSO_log2_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel log2FC_NVP2_vs_DMSO_minus -p max/2 -bs 25
computeMatrix scale-regions -S $log2FC_combo_vs_dTAG_plus_bw -R $SE_BED $other_BED -o $plus_reads_combo_vs_dTAG_log2_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel log2FC_NVP2_vs_DMSO_plus -p max/2 -bs 25
computeMatrix scale-regions -S $log2FC_combo_vs_dTAG_minus_bw -R $SE_BED $other_BED -o $minus_reads_combo_vs_dTAG_log2_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel log2FC_NVP2_vs_DMSO_minus -p max/2 -bs 25

## retain only signal from the respective gene's strand (and also antisense)
plus_strand_NVP2_vs_DMSO_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagenes_plus_strand_log2FC.matrix.gz
minus_strand_NVP2_vs_DMSO_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagenes_minus_strand_log2FC.matrix.gz
plus_antisense_NVP2_vs_DMSO_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagenes_plus_antisense_log2FC.matrix.gz
minus_antisense_NVP2_vs_DMSO_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagenes_minus_antisense_log2FC.matrix.gz
plus_strand_combo_vs_dTAG_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_plus_strand_log2FC.matrix.gz
minus_strand_combo_vs_dTAG_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_minus_strand_log2FC.matrix.gz
plus_antisense_combo_vs_dTAG_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_plus_antisense_log2FC.matrix.gz
minus_antisense_combo_vs_dTAG_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_minus_antisense_log2FC.matrix.gz

computeMatrixOperations filterStrand -m $plus_reads_NVP2_vs_DMSO_log2_matrix -o $plus_strand_NVP2_vs_DMSO_log2_matrix --strand +
computeMatrixOperations filterStrand -m $minus_reads_NVP2_vs_DMSO_log2_matrix -o $minus_strand_NVP2_vs_DMSO_log2_matrix --strand -
computeMatrixOperations filterStrand -m $plus_reads_NVP2_vs_DMSO_log2_matrix -o $plus_antisense_NVP2_vs_DMSO_log2_matrix --strand -
computeMatrixOperations filterStrand -m $minus_reads_NVP2_vs_DMSO_log2_matrix -o $minus_antisense_NVP2_vs_DMSO_log2_matrix --strand +
computeMatrixOperations filterStrand -m $plus_reads_combo_vs_dTAG_log2_matrix -o $plus_strand_combo_vs_dTAG_log2_matrix --strand +
computeMatrixOperations filterStrand -m $minus_reads_combo_vs_dTAG_log2_matrix -o $minus_strand_combo_vs_dTAG_log2_matrix --strand -
computeMatrixOperations filterStrand -m $plus_reads_combo_vs_dTAG_log2_matrix -o $plus_antisense_combo_vs_dTAG_log2_matrix --strand -
computeMatrixOperations filterStrand -m $minus_reads_combo_vs_dTAG_log2_matrix -o $minus_antisense_combo_vs_dTAG_log2_matrix --strand +

## merge the strands back together to get the "sense signal" matrix
merged_sense_NVP2_vs_DMSO_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagenes_merged_sense_log2FC.matrix.gz
merged_antisense_NVP2_vs_DMSO_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagenes_merged_antisense_log2FC.matrix.gz
merged_sense_combo_vs_dTAG_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_merged_sense_log2FC.matrix.gz
merged_antisense_combo_vs_dTAG_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_merged_antisense_log2FC.matrix.gz

computeMatrixOperations rbind -m $plus_strand_NVP2_vs_DMSO_log2_matrix $minus_strand_NVP2_vs_DMSO_log2_matrix -o $merged_sense_NVP2_vs_DMSO_log2_matrix
computeMatrixOperations rbind -m $plus_antisense_NVP2_vs_DMSO_log2_matrix $minus_antisense_NVP2_vs_DMSO_log2_matrix -o $merged_antisense_NVP2_vs_DMSO_log2_matrix
computeMatrixOperations rbind -m $plus_strand_combo_vs_dTAG_log2_matrix $minus_strand_combo_vs_dTAG_log2_matrix -o $merged_sense_combo_vs_dTAG_log2_matrix
computeMatrixOperations rbind -m $plus_antisense_combo_vs_dTAG_log2_matrix $minus_antisense_combo_vs_dTAG_log2_matrix -o $merged_antisense_combo_vs_dTAG_log2_matrix

## finally plot the sense PRO-seq signal heatmaps & metagenes
sense_NVP2_vs_DMSO_log2_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagenes_merged_sense_log2FC_heatmap.pdf
antisense_NVP2_vs_DMSO_log2_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagenes_merged_antisense_log2FC_heatmap.pdf
sense_combo_vs_dTAG_log2_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_merged_sense_log2FC_heatmap.pdf
antisense_combo_vs_dTAG_log2_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_merged_antisense_log2FC_heatmap.pdf

plotHeatmap -m $merged_sense_NVP2_vs_DMSO_log2_matrix -o $sense_NVP2_vs_DMSO_log2_output --dpi 300 --colorList "#050596,#ffffff00,#960505" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortRegions ascend #--zMax 1 --zMin "-1"
plotHeatmap -m $merged_antisense_NVP2_vs_DMSO_log2_matrix -o $antisense_NVP2_vs_DMSO_log2_output --dpi 300 --colorList "#050596,#ffffff00,#960505" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortRegions ascend #--zMax 1 --zMin "-1"
plotHeatmap -m $merged_sense_combo_vs_dTAG_log2_matrix -o $sense_combo_vs_dTAG_log2_output --dpi 300 --colorList "#050596,#ffffff00,#960505" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortRegions ascend #--zMax 1 --zMin "-1"
plotHeatmap -m $merged_antisense_combo_vs_dTAG_log2_matrix -o $antisense_combo_vs_dTAG_log2_output --dpi 300 --colorList "#050596,#ffffff00,#960505" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortRegions ascend #--zMax 1 --zMin "-1"


#############################################
##### cbind DMSO | 2h | log2FC heatmaps #####
#############################################
## merge left-to-right sample signal and log2FC

"""
check that the number of regions is the same:
zcat matrix.gz | wc -l
if yes, go ahead and cbind!
if no, remove regions from longer matrices using sort command.
"""
#9483 sense_combo_vs_dTAG_log2_matrix.bed       >>> this is the matrix with fewest regions
#9487 sense_NVP2_vs_log2_matrix.bed
#9520 sense_regions_4samples.bed
merged_sense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_merged_sense.matrix.gz
merged_sense_NVP2_vs_DMSO_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_NVP2-over-DMSO_metagenes_merged_sense_log2FC.matrix.gz
merged_sense_combo_vs_dTAG_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_merged_sense_log2FC.matrix.gz

# get the regions BED file for the matrix with fewest regions
short_matrix1_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_combo-over-dTAG_metagenes_merged_sense_log2FC.bed
plotHeatmap -m $merged_sense_combo_vs_dTAG_log2_matrix -o test2.pdf --sortRegions keep --outFileSortedRegions $short_matrix1_BED

# downsample the two longer matrices using the BED regions file from the shortest matrix
computeMatrixOperations sort -m $merged_sense_matrix -R $short_matrix1_BED -o temporary_clean_4samples.matrix.gz

# merge the first two matrices which now contain same regions
computeMatrixOperations cbind -m temporary_clean_4samples.matrix.gz $merged_sense_combo_vs_dTAG_log2_matrix -o temporary_clean_4samples_comboLog2FC.matrix.gz

# get regions from that 5-sample matrix
short_matrix2_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/temporary_clean_4samples_comboLog2FC.bed
plotHeatmap -m temporary_clean_4samples_comboLog2FC.matrix.gz -o test3.pdf --sortRegions keep --outFileSortedRegions $short_matrix2_BED

# downsample the last matrix using the BED regions file from the 5-sample matrix => some regions are not
computeMatrixOperations sort -m $merged_sense_NVP2_vs_DMSO_log2_matrix -R $short_matrix2_BED -o temporary_clean_NVP2_vs_DMSO_log2FC.matrix.gz

# and Vice versa (getting the regions from the downsampled 1-sample matrix first)
short_matrix3_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/temporary_clean_NVP2_vs_DMSO_log2FC.bed
plotHeatmap -m temporary_clean_NVP2_vs_DMSO_log2FC.matrix.gz -o test4.pdf --sortRegions keep --outFileSortedRegions $short_matrix3_BED
computeMatrixOperations sort -m temporary_clean_4samples_comboLog2FC.matrix.gz -R $short_matrix3_BED -o temporary_clean_5samples.matrix.gz

# merge the cleaned up 5-sample matrix and the 1-sample matrix which now contain same regions
sense_master_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_all-in-one_merged_sense_log2FC.matrix.gz
computeMatrixOperations cbind -m temporary_clean_5samples.matrix.gz temporary_clean_NVP2_vs_DMSO_log2FC.matrix.gz -o $sense_master_matrix

# rename and reorder the samples,
computeMatrixOperations relabel -m $sense_master_matrix -o $sense_master_matrix --sampleLabels "MJ-19-30_DMSO_2h" "MJ-19-30_NVP2_30min" "MJ-19-30_dTAG7_2h" "MJ-19-30_dTAG7_NVP2_combined" "log2FC_combo_vs_dTAG" "log2FC_NVP2_vs_DMSO"
computeMatrixOperations subset -m $sense_master_matrix --samples "MJ-19-30_DMSO_2h" "MJ-19-30_NVP2_30min" "log2FC_NVP2_vs_DMSO" "MJ-19-30_dTAG7_2h" "MJ-19-30_dTAG7_NVP2_combined" "log2FC_combo_vs_dTAG" -o $sense_master_matrix

## plot combined DMSO | 1h | log2FC signal heatmaps
sense_master_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_metagenes_all-in-one_merged_sense_log2FC_heatmap.pdf
plotHeatmap -m $sense_master_matrix -o $sense_master_output --dpi 300 --colorList "#ffffff00,#960505" "#ffffff00,#960505" "#0b3c7d,#231f20,#f2de05" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortUsingSamples 1 --zMin 0 0 "-0.5" --zMax 0.4 0.4 0.5



"""
missing:    now change the order of samples as desired and plot. Repeat same procedure for antisense matrix if desired.

#clean up temporary matrices
rm

"""


##############################################################
##### make TSS+/-500 matrices for NVP2 combo experiments #####
##############################################################
DMSO_2h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_DMSO_2h_merged_plus_dm6_normalized.bw
DMSO_2h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_DMSO_2h_merged_minus_dm6_normalized.bw
NVP2_30min_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/bigWig/MJ-19-30_NVP2_30min_merged_plus_dm6_normalized.bw
NVP2_30min_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/bigWig/MJ-19-30_NVP2_30min_merged_minus_dm6_normalized.bw
dTAG_2h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/bigWig/MJ-19-30_dTAG7_2h_merged_plus_dm6_normalized.bw
dTAG_2h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/bigWig/MJ-19-30_dTAG7_2h_merged_minus_dm6_normalized.bw
dTAG_NVP2_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/bigWig/MJ-19-30_dTAG7_NPV2_combo_merged_plus_dm6_normalized.bw
dTAG_NVP2_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/bigWig/MJ-19-30_dTAG7_NPV2_combo_merged_minus_dm6_normalized.bw

## computeMatrix for the plus and minus read signals separately
plus_reads_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_plus_reads_autoreg.matrix.gz
minus_reads_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_minus_reads_autoreg.matrix.gz

computeMatrix reference-point -S $DMSO_2h_plus_bw $NVP2_30min_plus_bw $dTAG_2h_plus_bw $dTAG_NVP2_plus_bw -R $SE_BED $autoreg_BED $other_BED -o $plus_reads_matrix -a 500 -b 500 --missingDataAsZero --skipZeros --samplesLabel MJ-19-30_DMSO_2h_plus MJ-19-30_NVP2_30min_plus MJ-19-30_dTAG7_2h_plus MJ-19-30_dTAG7_NVP2_combined_plus -p max/2
computeMatrix reference-point -S $DMSO_2h_minus_bw $NVP2_30min_minus_bw $dTAG_2h_minus_bw $dTAG_NVP2_minus_bw -R $SE_BED $autoreg_BED $other_BED -o $minus_reads_matrix -a 500 -b 500 --missingDataAsZero --skipZeros --samplesLabel MJ-19-30_DMSO_2h_minus MJ-19-30_NVP2_30min_minus MJ-19-30_dTAG7_2h_minus MJ-19-30_dTAG7_NVP2_combined_minus -p max/2 --scale "-1"

## retain only signal from the respective gene's strand (and also antisense)
plus_strand_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_plus_strand_autoreg.matrix.gz
minus_strand_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_minus_strand_autoreg.matrix.gz
plus_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_plus_antisense_autoreg.matrix.gz
minus_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_minus_antisense_autoreg.matrix.gz

computeMatrixOperations filterStrand -m $plus_reads_matrix -o $plus_strand_matrix --strand +
computeMatrixOperations filterStrand -m $minus_reads_matrix -o $minus_strand_matrix --strand -
computeMatrixOperations filterStrand -m $plus_reads_matrix -o $plus_antisense_matrix --strand -
computeMatrixOperations filterStrand -m $minus_reads_matrix -o $minus_antisense_matrix --strand +

## merge the strands back together to get the "sense signal" matrix
merged_sense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_sense_autoreg.matrix.gz  ### 9344 rows
merged_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_antisense_autoreg.matrix.gz ### 9243 rows

computeMatrixOperations rbind -m $plus_strand_matrix $minus_strand_matrix -o $merged_sense_matrix
computeMatrixOperations rbind -m $plus_antisense_matrix $minus_antisense_matrix -o $merged_antisense_matrix

## finally plot the sense PRO-seq signal heatmaps & metagenes
sense_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_sense_heatmap.pdf
antisense_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_antisense_heatmap.pdf

plotHeatmap -m $merged_sense_matrix -o $sense_output --dpi 300 --colorList "#ffffff00,#960505" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortUsingSamples 1 #--outFileSortedRegions sense_regions_4samples.bed --zMax 0.2 --zMin "-0.2"
plotHeatmap -m $merged_antisense_matrix -o $antisense_output --dpi 300 --colorList "#ffffff00,#050596" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortUsingSamples 1 #--zMax 0.2 --zMin "-0.2"

##################################################################################################
##### plot sense and anti merged left-to-right (need to clean up beforehand)
# get the regions BED file for the matrix with fewest regions
short_matrix1_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_antisense_regions.bed
plotHeatmap -m $merged_antisense_matrix -o test5.pdf --sortRegions keep --outFileSortedRegions $short_matrix1_BED

# downsample the two longer matrix using the BED regions file from the shorter matrix
computeMatrixOperations sort -m $merged_sense_matrix -R $short_matrix1_BED -o temporary_clean_sense.matrix.gz  ### this is now with 9116 regions

# get the regions BED file for the matrix with fewest regions
plotHeatmap -m temporary_clean_sense.matrix.gz -o test6.pdf --sortRegions keep --outFileSortedRegions temporary_clean_sense_regions.bed

# downsample the shorter matrix using the BED regions file from the cleaned-up longer matrix
computeMatrixOperations sort -m $merged_antisense_matrix -R temporary_clean_sense_regions.bed -o temporary_clean_antisense.matrix.gz  ### this is now with 9116 regions

# merge the two matrices which now contain same regions
merged_combined_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_both-senses.matrix.gz
computeMatrixOperations cbind -m temporary_clean_sense.matrix.gz temporary_clean_antisense.matrix.gz -o $merged_combined_matrix


## finally plot the sense PRO-seq signal heatmaps & metagenes
combined_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_both-senses_heatmap.pdf

plotHeatmap -m $merged_combined_matrix -o $combined_output --dpi 300 --colorList "#ffffff00,#960505" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortUsingSamples 1 #--outFileSortedRegions sense_regions_4samples.bed --zMax 0.2 --zMin "-0.2"
##################################################################################################

#######################################
##### MJ-19-30 TSS+/-500 profiles #####
#######################################
sense_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_sense_profile.pdf
antisense_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_antisense_profile.pdf

plotProfile -m $merged_sense_matrix -o $sense_profile_output --dpi 300 --plotWidth 10 --plotHeight 10 --colors "#be1e2d" "#231f20" --plotType se
plotProfile -m $merged_antisense_matrix -o $antisense_profile_output --dpi 300 --plotWidth 10 --plotHeight 10 --colors "#be1e2d" "#231f20" --plotType se

## now seperate by SE vs. non-SE rather than condition
sense_perGroup_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_sense_perGroup_profile_autoreg.pdf
antisense_perGroup_profile_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_TSS-region_merged_antisense_perGroup_profile.pdf

# direct plotting
plotProfile -m $merged_sense_matrix -o $sense_perGroup_profile_output --dpi 300 --plotWidth 10 --plotHeight 10 --colors "#231f2080" "#94919180" "#ed670780" "#be1e2d80" --plotType fill --perGroup
plotProfile -m $merged_antisense_matrix -o $antisense_perGroup_profile_output --dpi 300 --plotWidth 10 --plotHeight 10 --colors "#231f2080" "#94919180" "#ed670780" "#be1e2d80" --plotType fill --perGroup

# raw values for manual plotting
plotProfile -m $merged_sense_matrix -o $sense_perGroup_profile_output --dpi 300 --plotWidth 10 --plotHeight 10 --colors "#231f2080" "#94919180" "#ed670780" "#be1e2d80" --plotType fill --perGroup --outFileNameData MJ-19-30_TSS_autoregTF_6cell-smoothed_profile_values_raw_output.txt





















"""
## compute the master matrix, containing individual signal and log2FC
DMSO_2h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_DMSO_2h_merged_plus_dm6_normalized.bw
DMSO_2h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_DMSO_2h_merged_minus_dm6_normalized.bw
NVP2_30min_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/bigWig/MJ-19-30_NVP2_30min_merged_plus_dm6_normalized.bw
NVP2_30min_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/bigWig/MJ-19-30_NVP2_30min_merged_minus_dm6_normalized.bw
dTAG_2h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/bigWig/MJ-19-30_dTAG7_2h_merged_plus_dm6_normalized.bw
dTAG_2h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/bigWig/MJ-19-30_dTAG7_2h_merged_minus_dm6_normalized.bw
dTAG_NVP2_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/bigWig/MJ-19-30_dTAG7_NPV2_combo_merged_plus_dm6_normalized.bw
dTAG_NVP2_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/bigWig/MJ-19-30_dTAG7_NPV2_combo_merged_minus_dm6_normalized.bw
log2FC_NVP2_vs_DMSO_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_log2FC_NVP2-over-DMSO_plus_normalized.bw
log2FC_NVP2_vs_DMSO_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_log2FC_NVP2-over-DMSO_minus_normalized.bw
log2FC_combo_vs_dTAG_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_log2FC_combo_vs_dTAG_plus_normalized.bw
log2FC_combo_vs_dTAG_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bigWig/MJ-19-30_log2FC_combo_vs_dTAG_minus_normalized.bw

plus_reads_master_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_all-in-one_metagenes_plus_reads_log2FC.matrix.gz
minus_reads_master_log2_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_all-in-one_metagene_minus_reads_log2FC.matrix.gz

computeMatrix scale-regions -S $DMSO_2h_plus_bw $NVP2_30min_plus_bw $log2FC_NVP2_vs_DMSO_plus_bw $dTAG_2h_plus_bw $dTAG_NVP2_plus_bw $log2FC_combo_vs_dTAG_plus_bw -R $SE_BED $other_BED -o $plus_reads_master_log2_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel DMSO_2h NVP2_30min log2FC_NVP2_vs_DMSO dTAG7_2h dTAG7_2h_NVP2_30min log2FC_combo_vs_dTAG -p max/2 -bs 25
computeMatrix scale-regions -S $DMSO_2h_minus_bw $NVP2_30min_minus_bw $log2FC_NVP2_vs_DMSO_minus_bw $dTAG_2h_minus_bw $dTAG_NVP2_minus_bw $log2FC_combo_vs_dTAG_minus_bw -R $SE_BED $other_BED -o $minus_reads_master_log2_matrix -b 2000 -m 10000 -a 5000 --unscaled5prime 500 --unscaled3prime 500 --missingDataAsZero --skipZeros --samplesLabel DMSO_2h NVP2_30min log2FC_NVP2_vs_DMSO dTAG7_2h dTAG7_2h_NVP2_30min log2FC_combo_vs_dTAG -p max/2 -bs 25

## retain only signal from the respective gene's strand (and also antisense)
plus_strand_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_all-in-one_metagenes_plus_strand.matrix.gz
minus_strand_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_all-in-one_metagenes_minus_strand.matrix.gz
plus_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_all-in-one_metagenes_plus_antisense.matrix.gz
minus_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_all-in-one_metagenes_minus_antisense.matrix.gz

computeMatrixOperations filterStrand -m $plus_reads_master_log2_matrix -o $plus_strand_matrix --strand +
computeMatrixOperations filterStrand -m $minus_reads_master_log2_matrix -o $minus_strand_matrix --strand -
computeMatrixOperations filterStrand -m $plus_reads_master_log2_matrix -o $plus_antisense_matrix --strand -
computeMatrixOperations filterStrand -m $minus_reads_master_log2_matrix -o $minus_antisense_matrix --strand +

## merge the strands back together to get the "sense signal" matrix
merged_sense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_all-in-one_metagenes_merged_sense.matrix.gz
merged_antisense_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_all-in-one_metagenes_merged_antisense.matrix.gz

computeMatrixOperations rbind -m $plus_strand_matrix $minus_strand_matrix -o $merged_sense_matrix
computeMatrixOperations rbind -m $plus_antisense_matrix $minus_antisense_matrix -o $merged_antisense_matrix

## finally plot the sense PRO-seq signal heatmaps & metagenes
sense_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_all-in-one_metagenes_merged_sense_heatmap.pdf
antisense_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-30_all-in-one_metagenes_merged_antisense_heatmap.pdf

plotHeatmap -m $merged_sense_matrix -o $sense_output --dpi 300 --colorList "#ffffff00,#960505" "#ffffff00,#960505" "#e66101,#f7f7f7,#5e3c99" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortRegions keep --zMin 0 0 "-1" --zMax 0.002 0.002 1 --outFileSortedRegions sense_regions_all-in-one_manually_sorted.bed
plotHeatmap -m $merged_antisense_matrix -o $antisense_output --dpi 300 --colorList "#ffffff00,#960505" "#ffffff00,#960505" "#5e3c99,#f7f7f7,#e66101" --alpha 1 --regionsLabel "464 SEs" "9276 non-SEs" --sortUsing mean --sortUsingSamples 1 --zMin 0 0 "-1" --zMax 0.002 0.002 1
"""
