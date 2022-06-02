# Parse arguments for chipseq data path
while getopts 'd:' flag
do 
    case "${flag}" in 
        (d) chipseq_PATH=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (d) exit 9999;;
            esac;;
    esac
done

chipseq_SAMPLESHEET=$chipseq_PATH'/samplesheet.csv'
chipseq_OUTPUT=$chipseq_PATH'/chipseq_output';
echo '========== Initialized DHM workflow =========='

# Run nf-core/chip_seq
echo 'Create output folder'
mkdir -p $chipseq_OUTPUT;
echo 'Make sample sheet'
python utils/fastq_dir_to_samplesheet.py $chipseq_PATH $chipseq_SAMPLESHEET chipseq -r1 '_1.fastq.gz' -r2 '_2.fastq.gz' -ctl MCF7_DMSO -ab 2 -ri 3 -gi 0,1,2;
echo 'Run nf-core/chipseq'
nextflow run nf-core/chipseq --input $chipseq_SAMPLESHEET --outdir $chipseq_OUTPUT --genome GRCh38 -profile docker --narrow_peak true;
echo '========== Finished nf-core/chipseq =========='
