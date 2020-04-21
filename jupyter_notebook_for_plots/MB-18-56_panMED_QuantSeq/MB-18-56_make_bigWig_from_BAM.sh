#!/bin/bash
#SBATCH --output /scratch/lab_winter/martin/playground/test/MB-18-58_bam_to_bw.%j.log
#SBATCH --job-name=bam_to_bw
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
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
ref_index_sizes=/data/groups/lab_winter/reference_files/genomes/hg38_SIRV.chrom.sizes

######### index the aligned bam files
#find /scratch/lab_winter/matthias/MB18_56_QuantSeq_merge_analysis/results/results_pipeline/*/*.bam | parallel --no-notice 'samtools index {}'

######### run deeptools bamCoverage separately on strands
#find /scratch/lab_winter/matthias/MB18_56_QuantSeq_merge_analysis/results/results_pipeline/*/*.bam | parallel --no-notice 'bamCoverage -b {} -o {.}_plus_CPM.bw -p max/2 -bs 1 --normalizeUsing CPM --effectiveGenomeSize 2913022398 --samFlagExclude 16 -v'
#find /scratch/lab_winter/matthias/MB18_56_QuantSeq_merge_analysis/results/results_pipeline/*/*.bam | parallel --no-notice 'bamCoverage -b {} -o {.}_minus_CPM.bw -p max/2 -bs 1 --normalizeUsing CPM --effectiveGenomeSize 2913022398 --samFlagInclude 16 -v'

# bedgraph
find /scratch/lab_winter/matthias/MB18_56_QuantSeq_merge_analysis/results/results_pipeline/*/*.bam | parallel --no-notice "genomeCoverageBed -bg -strand '+' -ibam {} -g $ref_index_sizes | LC_COLLATE=C sort -k1,1 -k2,2n - > {.}_plus_unnorm.bg"
find /scratch/lab_winter/matthias/MB18_56_QuantSeq_merge_analysis/results/results_pipeline/*/*.bam | parallel --no-notice "genomeCoverageBed -bg -scale -1 -strand '-' -ibam {} -g $ref_index_sizes | LC_COLLATE=C sort -k1,1 -k2,2n - > {.}_minus_unnorm.bg"

# bigWig
find /scratch/lab_winter/matthias/MB18_56_QuantSeq_merge_analysis/results/results_pipeline/*/*_plus_unnorm.bg | parallel --no-notice "bedGraphToBigWig {} $ref_index_sizes {.}.bw"
find /scratch/lab_winter/matthias/MB18_56_QuantSeq_merge_analysis/results/results_pipeline/*/*_minus_unnorm.bg | parallel --no-notice "bedGraphToBigWig {} $ref_index_sizes {.}.bw"
