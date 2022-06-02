checkIfDirExist(){
    allPaths=("$@")
    for Path_ in "${allPaths[@]}"
    do
        if [ ! -d $Path_ ] 
            then
                echo 'Please check if input folders are named '$Path_
                exit 9999
            fi
    done
    }
    
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

mrnaseq_PATH=$INPUT_PATH'/mRNA_seq';
RNAseq_SAMPLESHEET=$mrnaseq_PATH'/samplesheet.csv'
RNAseq_OUTPUT=$mrnaseq_PATH'/rnaseq_output';

chipseq_PATH=$INPUT_PATH'/ChIP_seq';
chipseq_SAMPLESHEET=$chipseq_PATH'/samplesheet.csv'
chipseq_OUTPUT=$chipseq_PATH'/chipseq_output';
allPaths=( $mrnaseq_PATH $chipseq_PATH )
checkIfDirExist "${allPaths[@]}"

# Run nf-core/rna_seq
# echo 'Create output folder'
# mkdir -p $RNAseq_OUTPUT;
# echo 'Make sample sheet'
# python $nfcore_rnaseq/fastq_dir_to_samplesheet.py $mrnaseq_PATH $RNAseq_SAMPLESHEET rnaseq -st unstranded -r1 '_1.fastq.gz' -r2 '_2.fastq.gz';
# echo 'Run nf-core/rnaseq'
# nextflow run nf-core/rnaseq --input $RNAseq_SAMPLESHEET --outdir $RNAseq_OUTPUT --genome GRCh38 -profile docker
# echo 'Finished.'

# Run nf-core/chip_seq
echo 'Create output folder'
mkdir -p $chipseq_OUTPUT;
echo 'Make sample sheet'
python utils/fastq_dir_to_samplesheet.py $chipseq_PATH $chipseq_SAMPLESHEET chipseq -r1 '_1.fastq.gz' -r2 '_2.fastq.gz' -ctl MCF7_DMSO -ab 2 -ri 3 -gi 0,1,2;
echo 'Run nf-core/chipseq'
nextflow run nf-core/chipseq --input $chipseq_SAMPLESHEET --outdir $chipseq_OUTPUT --genome GRCh38 -profile docker --narrow_peak true;
echo 'Finished.'
