# Parse arguments for input path, output path
while getopts 'g:d:' flag
do 
    case "${flag}" in 
        (g) GENE=${OPTARG};;
        (d) DATAPATH=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (g) exit 1;;
                (d) exit 1;;
            esac;;
    esac
done

OUTPUT_PATH='/home/dhthutrang/Krebs/episplice-pipeline/utils/result/'$GENE
mkdir $OUTPUT_PATH 
grep $GENE -w /home/dhthutrang/Krebs/episplice-pipeline/utils/prepare_refgen/refgen_collapsed.gtf | grep -v aggregate_gene > $OUTPUT_PATH/FLANK.txt
grep $GENE -w $DATAPATH/mRNA_seq/DEXSeq_output/csv/*_res.csv > $OUTPUT_PATH/DEU.txt

for file in $DATAPATH/ChIP_seq/manorm_output_annotated/*.bed
do
    echo $file
    FILENAME=${file##*/}
    FILENAME=${FILENAME%%_*}
    grep $GENE -w $file > $OUTPUT_PATH/$FILENAME'_DHM.txt'

done
