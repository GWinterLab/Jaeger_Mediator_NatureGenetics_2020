################################################################
##### dREG was run on MJ-19-14_DMSO_1h_merged bigWig files #####
################################################################

SE_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/KBM7_H3K27ac_hg38_peaks_Gateway_SuperEnhancers.bed # wc -l is 532
dREG_peaks_full=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/dREG_results/MJ-19-14_DMSO_1h_merged.dREG.peak.full.bed # wc -l is 58868

dREG_filtered_centers=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/dREG_results/MJ-19-14_DMSO_1h_merged.dREG.peak.filtered.centers.bed

# filter dREG peaks for score > 0.5; output BED file with 1bp center positions (wc -l is 37790)
awk 'BEGIN{OFS="\t"} ($4 > 0.5) {print $1,$6,$6+1,"name",$4,"."}' $dREG_peaks_full > $dREG_filtered_centers


dREG_SE_constituents=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/dREG_results/MJ-19-14_DMSO_1h_merged.dREG.peak.filtered.centers.SEconstituents.bed

# intersect dREG centers that overlap SEs (wc -l is 2760)
intersectBed -a $dREG_filtered_centers -b $SE_BED -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$10,$5,$6}' > $dREG_SE_constituents


ROSE_input_enhancers=/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/KBM7_H3K27ac_hg38_peaks.noOverlap_H3K4me3_TSS.bed

dREG_nonSE_enhancers=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/dREG_results/MJ-19-14_DMSO_1h_merged.dREG.peak.filtered.centers.nonSEenhancers.bed

# intersect non-SE dREG enhancers (filtered H3K27ac peaks) (wc -l is 6160)
intersectBed -a $dREG_filtered_centers -b $SE_BED -v | intersectBed -a stdin -b $ROSE_input_enhancers -wb | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$10,$5,$6}' > $dREG_nonSE_enhancers


#################################################################
##### use deeptools 3.3.0 to make the matrices and heatmaps #####
#################################################################
DMSO_1h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/bigWig/MJ-19-14_DMSO_1h_merged_plus_dm6_normalized.bw
DMSO_1h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/bigWig/MJ-19-14_DMSO_1h_merged_minus_dm6_normalized.bw
dTAG_1h_plus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/bigWig/MJ-19-14_dTAG7_1h_merged_plus_dm6_normalized.bw
dTAG_1h_minus_bw=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/bigWig/MJ-19-14_dTAG7_1h_merged_minus_dm6_normalized.bw
K27ac_bw=/scratch/lab_winter/hana/2_Mediator_SE_CRC/KBM7_cell_line/mapping_hg38/KBM7_H3K27ac_hg38.bigWig

combined_matrix=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_combined_K27ac.matrix.gz

combined_output=/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/enrichment_heatmaps/MJ-19-14_1h_combined_matrix_K27ac.pdf

computeMatrix reference-point -S $DMSO_1h_plus_bw $dTAG_1h_plus_bw $DMSO_1h_minus_bw $dTAG_1h_minus_bw $K27ac_bw -R $dREG_SE_constituents $dREG_nonSE_enhancers -o $combined_matrix -a 1000 -b 1000 --missingDataAsZero --skipZeros --samplesLabel DMSO_1h_plus dTAG7_1h_plus DMSO_1h_minus dTAG7_1h_minus KBM7_H3K27ac




plotHeatmap -m $combined_matrix -o $combined_output --dpi 300 --colorList "#050596,#ffffff00,#960505" --alpha 1 --zMax 0.2 --zMin "-0.2" --refPointLabel "center" --regionsLabel "SE constituents" "other dREG enhancers" --sortUsing max --sortUsingSamples 1
