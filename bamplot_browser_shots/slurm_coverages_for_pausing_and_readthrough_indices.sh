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

#########################################
##### make all the BED region files #####
#########################################

NM_TT_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/max_RPKM_per_gene.coverage-TTseq_KBM7_DMSO1h2h_read1_merged.extended100bp.RPKM.only_NM.bed
NM_full_length=/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/BED_TT_expressed_genes/9739_full_length_expressed.bed
NM_TSS_window=/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/BED_TT_expressed_genes/9438_TSS-50_to_TSS+250_TSSwindow_expressed.bed
NM_genebody=/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/BED_TT_expressed_genes/9438_TSS+500_to_TTS_genebody_expressed.bed
NM_gene_end=/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/BED_TT_expressed_genes/9438_TTS-1000_to_TTS_gene_end_expressed.bed
NM_termination_window=/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/BED_TT_expressed_genes/9438_TTS+3000_to_TTS+6000_termination_window_expressed.bed

awk '($13 > 0.110155)' $NM_TT_cov | cut -f1-6 | sort -k1,1 -k2,2n > $NM_full_length #cutoff = median TT-seq RPKM

# define new bed window from -50 to +250 of the TSS (only genes >= 2000 bp)
awk 'BEGIN{OFS="\t"} ($6 == "+") && ($3-$2 >= 2000) {print $1,$2-50,$2+250,$4,$5,$6}; ($6 == "-") && ($3-$2 >= 2000) {print $1,$3-250,$3+50,$4,$5,$6}' $NM_full_length | sort -k1,1 -k2,2n > $NM_TSS_window
# define new genebody window from TSS+500 to the gene end (only genes >= 2000 bp)
awk 'BEGIN{OFS="\t"} ($6 == "+") && ($3-$2 >= 2000) {print $1,$2+500,$3,$4,$5,$6}; ($6 == "-") && ($3-$2 >= 2000) {print $1,$2,$3-500,$4,$5,$6}' $NM_full_length | sort -k1,1 -k2,2n > $NM_genebody

# define new termination window from TTS+3000 to TTS+6000 (only genes >= 2000 bp and not overlapping any other gene)
awk 'BEGIN{OFS="\t"} ($6 == "+") && ($3-$2 >= 2000) {print $1,$3+3000,$3+6000,$4,$5,$6}; ($6 == "-") && ($3-$2 >= 2000) {print $1,$2-6000,$2-3000,$4,$5,$6}' $NM_full_length | sort -k1,1 -k2,2n > $NM_termination_window
intersectBed -a $NM_termination_window -b $NM_full_length -v > tmp.file
mv tmp.file $NM_termination_window
# define new gene end window from TTS-1000 to the gene end (only genes >= 2000 bp)
awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2-4000,$3-6000,$4,$5,$6}; ($6 == "-") {print $1,$2+6000,$3+4000,$4,$5,$6}' $NM_termination_window | sort -k1,1 -k2,2n > $NM_gene_end


###########################################
##### compute norm MJ-19-14 coverages #####
###########################################
DMSO_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/bed_1bp/MJ-19-14_DMSO_1h_merged_1bp.sorted.bed.gz
dTAG_1h_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/bed_1bp/MJ-19-14_dTAG7_1h_merged_1bp.sorted.bed.gz


# TSS_window
DMSO_TSS_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/refGene_coverage/MJ-19-14_DMSO_1h_merged_TSSwindow_expressed_norm.cov
dTAG_1h_TSS_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/refGene_coverage/MJ-19-14_dTAG7_1h_merged_TSSwindow_expressed_norm.cov
coverageBed -a $NM_TSS_window -b $DMSO_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.247928; print}' - > $DMSO_TSS_cov
coverageBed -a $NM_TSS_window -b $dTAG_1h_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.176419; print}' - > $dTAG_1h_TSS_cov

# TSS_genebody
DMSO_genebody_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/refGene_coverage/MJ-19-14_DMSO_1h_merged_genebody_expressed_norm.cov
dTAG_1h_genebody_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/refGene_coverage/MJ-19-14_dTAG7_1h_merged_genebody_expressed_norm.cov
coverageBed -a $NM_genebody -b $DMSO_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.247928; print}' - > $DMSO_genebody_cov
coverageBed -a $NM_genebody -b $dTAG_1h_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.176419; print}' - > $dTAG_1h_genebody_cov

# TTS_gene_end
DMSO_gene_end_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/refGene_coverage/MJ-19-14_DMSO_1h_merged_gene_end_expressed_norm.cov
dTAG_1h_gene_end_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/refGene_coverage/MJ-19-14_dTAG7_1h_merged_gene_end_expressed_norm.cov
coverageBed -a $NM_gene_end -b $DMSO_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.247928; print}' - > $DMSO_gene_end_cov
coverageBed -a $NM_gene_end -b $dTAG_1h_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.176419; print}' - > $dTAG_1h_gene_end_cov

# TTS_termination_window
DMSO_termination_window_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/refGene_coverage/MJ-19-14_DMSO_1h_merged_termination_window_expressed_norm.cov
dTAG_1h_termination_window_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/refGene_coverage/MJ-19-14_dTAG7_1h_merged_termination_window_expressed_norm.cov
coverageBed -a $NM_termination_window -b $DMSO_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.247928; print}' - > $DMSO_termination_window_cov
coverageBed -a $NM_termination_window -b $dTAG_1h_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.176419; print}' - > $dTAG_1h_termination_window_cov



###########################################
##### compute norm MJ-19-30 coverages #####
###########################################
DMSO_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/bed_1bp/MJ-19-30_DMSO_2h_merged_1bp.sorted.bed.gz
dTAG_2h_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/bed_1bp/MJ-19-30_dTAG7_2h_merged_1bp.sorted.bed.gz
NVP2_30min_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/bed_1bp/MJ-19-30_NVP2_30min_merged_1bp.sorted.bed.gz
combo_BED=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/bed_1bp/MJ-19-30_dTAG7_NPV2_combo_merged_1bp.sorted.bed.gz


# TSS_window
DMSO_TSS_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/refGene_coverage/MJ-19-30_DMSO_2h_merged_TSSwindow_expressed_norm.cov
dTAG_2h_TSS_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/refGene_coverage/MJ-19-30_dTAG7_2h_merged_TSSwindow_expressed_norm.cov
NVP2_30min_TSS_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/refGene_coverage/MJ-19-30_NVP2_30min_merged_TSSwindow_expressed_norm.cov
combo_TSS_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/refGene_coverage/MJ-19-30_dTAG7_NPV2_combo_merged_TSSwindow_expressed_norm.cov
coverageBed -a $NM_TSS_window -b $DMSO_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.412471; print}' - > $DMSO_TSS_cov
coverageBed -a $NM_TSS_window -b $dTAG_2h_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.274846; print}' - > $dTAG_2h_TSS_cov
coverageBed -a $NM_TSS_window -b $NVP2_30min_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.392713; print}' - > $NVP2_30min_TSS_cov
coverageBed -a $NM_TSS_window -b $combo_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.37138; print}' - > $combo_TSS_cov

# TSS_genebody
DMSO_genebody_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/refGene_coverage/MJ-19-30_DMSO_2h_merged_genebody_expressed_norm.cov
dTAG_2h_genebody_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/refGene_coverage/MJ-19-30_dTAG7_2h_merged_genebody_expressed_norm.cov
NVP2_30min_genebody_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/refGene_coverage/MJ-19-30_NVP2_30min_merged_genebody_expressed_norm.cov
combo_genebody_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/refGene_coverage/MJ-19-30_dTAG7_NPV2_combo_merged_genebody_expressed_norm.cov
coverageBed -a $NM_genebody -b $DMSO_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.412471; print}' - > $DMSO_genebody_cov
coverageBed -a $NM_genebody -b $dTAG_2h_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.274846; print}' - > $dTAG_2h_genebody_cov
coverageBed -a $NM_genebody -b $NVP2_30min_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.392713; print}' - > $NVP2_30min_genebody_cov
coverageBed -a $NM_genebody -b $combo_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.37138; print}' - > $combo_genebody_cov

# TTS_gene_end
DMSO_gene_end_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/refGene_coverage/MJ-19-30_DMSO_2h_merged_gene_end_expressed_norm.cov
dTAG_2h_gene_end_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/refGene_coverage/MJ-19-30_dTAG7_2h_merged_gene_end_expressed_norm.cov
NVP2_30min_gene_end_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/refGene_coverage/MJ-19-30_NVP2_30min_merged_gene_end_expressed_norm.cov
combo_gene_end_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/refGene_coverage/MJ-19-30_dTAG7_NPV2_combo_merged_gene_end_expressed_norm.cov
coverageBed -a $NM_gene_end -b $DMSO_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.412471; print}' - > $DMSO_gene_end_cov
coverageBed -a $NM_gene_end -b $dTAG_2h_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.274846; print}' - > $dTAG_2h_gene_end_cov
coverageBed -a $NM_gene_end -b $NVP2_30min_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.392713; print}' - > $NVP2_30min_gene_end_cov
coverageBed -a $NM_gene_end -b $combo_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.37138; print}' - > $combo_gene_end_cov

# TTS_termination_window
DMSO_termination_window_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/refGene_coverage/MJ-19-30_DMSO_2h_merged_termination_window_expressed_norm.cov
dTAG_2h_termination_window_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/refGene_coverage/MJ-19-30_dTAG7_2h_merged_termination_window_expressed_norm.cov
NVP2_30min_termination_window_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/refGene_coverage/MJ-19-30_NVP2_30min_merged_termination_window_expressed_norm.cov
combo_termination_window_cov=/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/refGene_coverage/MJ-19-30_dTAG7_NPV2_combo_merged_termination_window_expressed_norm.cov
coverageBed -a $NM_termination_window -b $DMSO_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.412471; print}' - > $DMSO_termination_window_cov
coverageBed -a $NM_termination_window -b $dTAG_2h_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.274846; print}' - > $dTAG_2h_termination_window_cov
coverageBed -a $NM_termination_window -b $NVP2_30min_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.392713; print}' - > $NVP2_30min_termination_window_cov
coverageBed -a $NM_termination_window -b $combo_BED -s -sorted | awk 'BEGIN{OFS="\t"} {$7 = $7 * 0.37138; print}' - > $combo_termination_window_cov
