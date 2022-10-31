# dir is the path to fastq_files folder
# metadata is the path to metadata.tsv file. Modify the last column of metadata.tsv to the new name
while getopts 'd:m:' flag
do 
    case "${flag}" in 
        (d) dir=${OPTARG};;
        (m) metadata=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (d) exit 1;;
                (m) exit 1;;
            esac;;
    esac
done
cd $dir/mRNA_seq/rnaseq_output/star_salmon
for file in *.bam
do
    echo $file 
    sample=${file%.markdup*}
    mv $file ${dir##*/}_${sample//_/}.bam
done
for file in *.bam.bai
do
    echo $file 
    sample=${file%.markdup*}
    mv $file ${dir##*/}_${sample//_/}.bam.bai
done