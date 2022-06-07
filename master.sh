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

rnaseq_PATH=$INPUT_PATH'/mRNA_seq';
chipseq_PATH=$INPUT_PATH'/ChIP_seq';
allPaths=( $mrnaseq_PATH $chipseq_PATH )
checkIfDirExist "${allPaths[@]}"

./DEU_workflow.sh -d $rnaseq_PATH
#./DHM_workflow.sh -d $chipseq_PATH