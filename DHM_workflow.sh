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

fastq_PATH=$chipseq_PATH'/fastq_files'
SAMPLESHEET=$chipseq_PATH'/samplesheet.csv'
chipseq_OUTPUT=$chipseq_PATH'/chipseq_output';
echo '========== Initialized DHM workflow =========='

# Run nf-core/chip_seq
echo 'Create output folder'
mkdir -p $chipseq_OUTPUT;
# echo 'Make sample sheet'
# python utils/fastq_dir_to_samplesheet.py $fastq_PATH $SAMPLESHEET chipseq -r1 '_1.fastq.gz' -r2 '_2.fastq.gz' -ctl MCF7_DMSO -ab 2 -ri 3 -gi 0,1,2; #Reddy
# python utils/fastq_dir_to_samplesheet.py $fastq_PATH $SAMPLESHEET chipseq -r1 '.fastq.gz' -r2 '.fastq.gz' -ctl MCF7_parental -ab 2 -ri 3 -gi 0,1,2; #Shen
# python utils/fastq_dir_to_samplesheet.py $fastq_PATH $SAMPLESHEET chipseq -r1 '.fastq.gz' -r2 '.fastq.gz' -ctl MCF7_ctrl -ab 2 -ri 3 -gi 0,1,2; #Shen
# echo 'Run nf-core/chipseq'
# nextflow run nf-core/chipseq --input $SAMPLESHEET --outdir $chipseq_OUTPUT --genome GRCh38 --save_reference true -profile docker --narrow_peak true; #Reddy 16.06.22
# nextflow run nf-core/chipseq --input $SAMPLESHEET --outdir $chipseq_OUTPUT --genome GRCh38 --save_reference true -profile docker --narrow_peak true --skip_fastqc true --skip_picard_metrics true --skip_preseq true --skip_plot_profile true --skip_plot_fingerprint true --skip_spp true --skip_igv true --skip_diff_analysis true --keep_dups true; #Reddy 17.06.22
# nextflow run nf-core/chipseq --input $SAMPLESHEET --outdir $chipseq_OUTPUT --genome GRCh38 --save_reference true -profile docker --narrow_peak true --single_end --skip_plot_fingerprint true --skip_plot_profile true --skip_picard_metrics true --skip_preseq true; #Shen
# nextflow run nf-core/chipseq --input $SAMPLESHEET --outdir $chipseq_OUTPUT --genome GRCh38 --save_reference true -profile docker --narrow_peak true --single_end; #Messier
echo '========== Finished nf-core/chipseq =========='
