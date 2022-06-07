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

fastq_PATH=$rnaseq_PATH/fastq_files;
rnaseq_OUTPUT=$rnaseq_PATH'/rnaseq_output';
SAMPLESHEET=$fastq_PATH'/samplesheet.csv';
SAM_PATH=$rnaseq_PATH/'SAM_files';
count_PATH=$rnaseq_PATH/'count_files';
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

echo 'Generate exon count'
# mkdir -p $count_PATH
DEU_scripts/generate_exon_count.sh -i $SAM_PATH -o $count_PATH -g $ENCODE_REFGEN/reference_genome.2021_.corrected.gtf

echo '========== Generated exon counts =========='