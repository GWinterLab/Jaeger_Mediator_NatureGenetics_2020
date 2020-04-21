########################
### MJ-19-30 PRO-seq ###
########################

#title
TITLE='MJ-19-30_TSS_FAM91A1_plus_1000bins'
#output directory
OUTPUT='/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/bamplot/MJ-19-30_PROseq/'
#full-length region
#CHR='chr8:+:123767000-123825000'
#TSS region
CHR='chr8:+:123767000-123771000'
#readthrough region
#CHR='chr8:+:123805000-123835000'
#sense
SENSE='+'


#input bams
PROBAM1='/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_DMSO_2h_merged/aligned/MJ-19-30_DMSO_2h_merged_1bp.sorted.bam'
PROBAM2='/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_2h_merged/aligned/MJ-19-30_dTAG7_2h_merged_1bp.sorted.bam'
PROBAM3='/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_dTAG7_NPV2_combo_merged/aligned/MJ-19-30_dTAG7_NPV2_combo_merged_1bp.sorted.bam'
PROBAM4='/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-30_NVP2_30min_merged/aligned/MJ-19-30_NVP2_30min_merged_1bp.sorted.bam'

NAME1='MJ-19-30_PRO-seq_DMSO_2h_merged'
NAME2='MJ-19-30_PRO-seq_dTAG7_2h_merged'
NAME3='MJ-19-30_PRO-seq_dTAG7_NVP2_combo'
NAME4='MJ-19-30_PRO-seq_NVP2_30min_merged'

SCALE_FACTORS='0.412471,0.274846,0.37138,0.392713'

#### COMMAND
if [ $SENSE = '+' ]
then
   COLOR='236,37,38:236,37,38:236,37,38:236,37,38'
else
   COLOR='47,53,142:47,53,142:47,53,142:47,53,142'
fi

bamplot -b $PROBAM1,$PROBAM2,$PROBAM3,$PROBAM4 -i $CHR -g hg38 -t $TITLE -y UNIFORM -o $OUTPUT -n $NAME1,$NAME2,$NAME3,$NAME4 -p multiple -c $COLOR -s $SENSE --scale $SCALE_FACTORS -r
