# bampath is the path to bamdirs folder
while getopts 'b:p:' flag
do 
    case "${flag}" in 
        (b) bamdirs=${OPTARG};;
        (p) prefix=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (b) exit 1;;
                (p) exit 1;;
            esac;;
    esac
done
cd $bamdirs
for file in *.bam
do
    echo $file 
    sample=${file%.markdup*}
    mv $file ${prefix}_${sample//_/}.bam
done
for file in *.bam.bai
do
    echo $file 
    sample=${file%.markdup*}
    mv $file ${prefix}_${sample//_/}.bam.bai
done

for file in *.majiq *.sj
do
    new_name=$(echo $file | sed 's/\.//') 
    mv $file $new_name
done
