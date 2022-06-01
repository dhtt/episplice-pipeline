# Parse arguments for input path, output path
while getopts 'i:' flag
do 
    case "${flag}" in 
        (i) INPUT_PATH=${OPTARG};;
        # (c) NEXTFLOW_CONFIG=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (i) exit 1;;
                # (c) exit 1;;
            esac;;
    esac
done

mRNA_seq_PATH=$INPUT_PATH'/mRNA_seq';
RNAseq_SAMPLESHEET=$mRNA_seq_PATH'/samplesheet.csv'
RNAseq_OUTPUT=$mRNA_seq_PATH'/rnaseq_output';

ChIP_seq_PATH=$INPUT_PATH'/ChIP_seq';
# TODO: Error if the folders containing fastq are not mRNA_seq and ChIP_seq

# $nfcore_rnaseq='/home/dhthutrang/tools/nf-core/rnaseq/bin'
echo 'Create output folder'
mkdir -p $RNAseq_OUTPUT;
echo 'Make sample sheet'
python $nfcore_rnaseq/fastq_dir_to_samplesheet.py $mRNA_seq_PATH $RNAseq_SAMPLESHEET -st unstranded -r1 '_1.fastq.gz' -r2 '_2.fastq.gz';
echo 'Run nf-core/rna_seq'
nextflow run nf-core/rnaseq --input $RNAseq_SAMPLESHEET --outdir $RNAseq_OUTPUT --genome GRCh38 -profile docker
echo 'Finished.'