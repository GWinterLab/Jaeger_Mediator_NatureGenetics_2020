#!/bin/bash
#SBATCH --output /scratch/lab_winter/martin/MED14_PRO-seq_final/code/slurm_deeptools_log2fc.%j.log
#SBATCH --job-name=deeptools
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=120000
#SBATCH --mail-type=end
#SBATCH --mail-user=mjaeger@cemm.oeaw.ac.at


echo "Enviromental variables"
echo "======================"

echo $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NAME
echo $SLURM_JOB_PARTITION
echo $SLURM_NTASKS
echo $SLURM_NPROCS
echo $SLURM_JOB_ID
echo $SLURM_JOB_NUM_NODES
echo $SLURM_NODELIST
echo $SLURM_CPUS_ON_NODE

echo "======================"

source /home/mjaeger/.bashrc

cd /home/mjaeger/

#################################################################
##### use deeptools 3.3.0 to make log2FC dTAG/DMSO heatmaps #####
#################################################################
SE_BED=/scratch/lab_winter/martin/Mediator_CRC/BED_for_deeptools/SEs_464_TTseq_NMonly.bed
other_BED=/scratch/lab_winter/martin/Mediator_CRC/BED_for_deeptools/not_SEs_9276_TTseq_NMonly.bed


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
