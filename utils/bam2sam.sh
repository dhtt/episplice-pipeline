# Parse arguments for input path, output path
while getopts 'i:o:' flag
do 
    case "${flag}" in 
        (i) INPUT_PATH=${OPTARG};;
        (o) OUTPUT_PATH=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (i) exit 9999;;
                (o) exit 9999;;
            esac;;
    esac
done

echo "bam2sam.sh: STARTED"
bam_to_sam(){
    ACC_NO=${f%%.*}
    ACC_NO=${ACC_NO##*/}
    echo "Converting $ACC_NO"
    (samtools sort -n -O SAM $f > $OUTPUT_PATH/$ACC_NO.sam) || echo "Err: Cannot convert BAM to SAM $ACC_NO"
}

for f in $INPUT_PATH/*.bam
do
    bam_to_sam $f &
done
wait

echo "bam2sam.sh: DONE"
