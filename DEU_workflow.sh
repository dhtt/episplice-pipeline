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

rnaseq_SAMPLESHEET=$rnaseq_PATH'/samplesheet.csv'
rnaseq_OUTPUT=$rnaseq_PATH'/rnaseq_output';
echo '========== Initialized DEU workflow =========='

# # Run nf-core/rna_seq
# echo 'Create output folder'
# mkdir -p $RNAseq_OUTPUT;
# echo 'Make sample sheet'
# python $nfcore_rnaseq/fastq_dir_to_samplesheet.py $mrnaseq_PATH $RNAseq_SAMPLESHEET rnaseq -st unstranded -r1 '_1.fastq.gz' -r2 '_2.fastq.gz';
# echo 'Run nf-core/rnaseq'
# nextflow run nf-core/rnaseq --input $RNAseq_SAMPLESHEET --outdir $RNAseq_OUTPUT --genome GRCh38 -profile docker
# echo '========== Finished nf-core/rnaseq =========='

# # Generate exon counts
# echo 'Convert BAM to SAM'
# mkdir -p $rnaseq_PATH/'SAM_files'
# utils/bam2sam.sh -i $rnaseq_OUTPUT'/star_salmon' -o $rnaseq_PATH/SAM_files

# echo 'Generate exon count'
# mkdir -p $rnaseq_PATH/'count_files'
# DEU_scripts/generate_exon_count.sh -i $rnaseq_PATH/SAM_files -o $rnaseq_PATH/count_files -g $ENCODE_REFGEN/reference_genome.gtf

echo '========== Generated exon counts =========='