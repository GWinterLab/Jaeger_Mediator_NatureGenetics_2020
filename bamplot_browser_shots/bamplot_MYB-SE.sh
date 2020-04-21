#title
TITLE='MJ-19-14_1h_MYB-SE_plus_1000bins'
#output directory
OUTPUT='/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/bamplot/1h_PROseq/'
#target region
CHR='chr6:+:135240000-135348000'
#sense
SENSE='+'


#input bams
PROBAM1='/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_DMSO_1h_merged/aligned/MJ-19-14_DMSO_1h_merged_1bp.sorted.bam'
PROBAM2='/scratch/lab_winter/martin/MED14_PRO-seq_final/results/MJ-19-14_dTAG7_1h_merged/aligned/MJ-19-14_dTAG7_1h_merged_1bp.sorted.bam'

NAME1='PRO-seq_MED14_DMSO_1h_merged'
NAME2='PRO-seq_MED14_dTAG7_1h_merged'

SCALE_FACTORS='0.247928,0.176419'

#### COMMAND
if [ $SENSE = '+' ]
then
   COLOR='236,37,38:236,37,38'
else
   COLOR='47,53,142:47,53,142'
fi

bamplot -b $PROBAM1,$PROBAM2 -i $CHR -g hg38 -t $TITLE -y UNIFORM -o $OUTPUT -n $NAME1,$NAME2 -p multiple -c $COLOR -s $SENSE --scale $SCALE_FACTORS -r


################
### ChIP-seq ###
################

K27ac_BAM='/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/KBM7_histone_ChIP-seq/KBM7_H3K27ac.trimmed.bowtie2.hg38.filtered.bam'
K4me3_BAM='/scratch/lab_winter/martin/MED14_PRO-seq_final/other_final_files/KBM7_histone_ChIP-seq/KBM7_H3K4me3.trimmed.bowtie2.hg38.filtered.bam'

#Edit the names to give each bam a title
NAME1='KBM7_H3K27ac'
NAME2='KBM7_H3K4me3'

#edit these variables to specify region, genome, title etc...
CHR='chr6:+:135240000-135348000'
GENOME='hg38'
OUTPUT='/scratch/lab_winter/martin/MED14_PRO-seq_final/plots/bamplot/1h_PROseq/'
TITLE='H3_ChIP-seq_MYB-SE'
YAXIS='UNIFORM'    #use either UNIFORM or RELATIVE
SENSE='both' #sense of reads plotted. use either '+', '-', or 'both'
PLOT='multiple' #use either 'multiple' or 'single'
COLOR='35,31,32:128,130,132'

bamplot -b $K27ac_BAM,$K4me3_BAM -i $CHR -g $GENOME -t $TITLE -y $YAXIS -o $OUTPUT -n $NAME1,$NAME2 -p $PLOT -c $COLOR -s $SENSE -r
