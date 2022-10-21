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
cd $dir
while IFS="," read -r rec1 rec2
do
    mv $rec1 $rec2
done < <(cut -f11,12 $metadata)