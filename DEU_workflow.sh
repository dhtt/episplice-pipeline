# Parse arguments for rnaseq data path
while getopts 'd:' flag
do 
    case "${flag}" in 
        (d) rnaseq_PATH=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (d) exit 9999;;
            esac;;
    esac
done

fastq_PATH=$rnaseq_PATH/fastq_files; #Files in fastq_PATH must follow naming format celltype_condition_repx_read.fastq.gz
rnaseq_OUTPUT=$rnaseq_PATH'/rnaseq_output';
refgen_PATH=$ENCODE_REFGEN/reference_genome.2021_.corrected.gtf
SAMPLESHEET=$fastq_PATH'/samplesheet.csv';
SAM_PATH=$rnaseq_PATH/'SAM_files'; # File name in count_PATH must be celltype_condition_repx.txt
count_PATH=$rnaseq_PATH/'count_files'; # File name in count_PATH must be celltype_condition_repx.txt
DEXSeq_output_PATH=$rnaseq_PATH'/DEXSeq_output'
echo '========== Initialized DEU workflow =========='

# # Run nf-core/rna_seq
# echo 'Create output folder'
# mkdir -p $rnaseq_OUTPUT;
# echo 'Make sample sheet'
# python $nfcore_rnaseq/fastq_dir_to_samplesheet.py $fastq_PATH $SAMPLESHEET rnaseq -st unstranded -r1 '_1.fastq.gz' -r2 '_2.fastq.gz';
# echo 'Run nf-core/rnaseq'
# nextflow run nf-core/rnaseq --input $SAMPLESHEET --outdir $rnaseq_OUTPUT --genome GRCh38 -profile docker
# echo '========== Finished nf-core/rnaseq =========='

# # Generate exon counts
# echo 'Convert BAM to SAM'
# mkdir -p $SAM_PATH
# utils/bam2sam.sh -i $rnaseq_OUTPUT'/star_salmon' -o $SAM_PATH

# echo 'Generate exon count'
# mkdir -p $count_PATH
# DEU_scripts/generate_exon_count.sh -i $SAM_PATH -o $count_PATH -g $refgen_PATH
# echo '========== Generated exon counts =========='

# DEXSeq run
echo 'Start DEXseq analysis'
mkdir -p $DEXSeq_output_PATH $DEXSeq_output_PATH/csv $DEXSeq_output_PATH/html $DEXSeq_output_PATH/r_data $DEXSeq_output_PATH/plot
Rscript DEU_scripts/DEXSeq_analysis.R -f $count_PATH -a MCF7_DMSO -b MCF7_50nM -G $refgen_PATH -n 8
echo '========== Generated exon counts =========='