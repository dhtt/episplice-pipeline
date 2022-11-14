# bampath is the path to bamdirs folder
while getopts 'b:' flag
do 
    case "${flag}" in 
        (b) bampath=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (b) exit 1;;
            esac;;
    esac
done
cd $bampath/bamdirs
for file in *.bam
do
    echo $file 
    sample=${file%.markdup*}
    mv $file ${bampath##*/}_${sample//_/}.bam
done
for file in *.bam.bai
do
    echo $file 
    sample=${file%.markdup*}
    mv $file ${bampath##*/}_${sample//_/}.bam.bai
done