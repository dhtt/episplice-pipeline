#FIRST ACTIVATE PYTHON 3.9 ENV source /home/dhthutrang/python3.9/bin/activate
while getopts 'g:p:' flag
do 
    case "${flag}" in 
        (g) groupid=${OPTARG};;
        (p) exppath=${OPTARG};;
        (:) 
            case ${OPTARG} in 
                (g) exit 1;;
                (p) exit 1;;
            esac;;
    esac
done


# Set up paths and dirs 
cd $exppath
group_ids=$groupid
sample_path=$exppath/$group_ids
bamdirs=$sample_path/bamdirs
refgen_path=/home/dhthutrang/ENCODE/refgen/majiq/hg38.ncbiRefSeq.gtf
build_result=$sample_path/build_result

echo 'PATH TO SAMPLES: '$sample_path
mkdir -p $bamdirs
IFS='_' read -r -a array <<< "$group_ids"


# Create a config file
config_path=$sample_path/config.ini
echo 'CREATE CONFIG FILE: '$config_path
printf "[info]\nbamdirs=$bamdirs\ngenome=hg38\n\n[experiments]\n" > $config_path
for index in "${!array[@]}"
do
    treatment=${array[index]}
    sample_ids=$bamdirs/$treatment*.bam
    printf "$treatment=" >> $config_path
    echo $sample_ids | sed -r 's#'"$bamdirs/"'##g' | sed -r 's#.bam #,#g' | sed -r 's#.bam##g' >> $config_path
done


# Build MAJIQ index
echo 'START MAJIQ BUILD'
majiq build --disable-ir --disable-denovo --disable-denovo-ir $refgen_path -j 8 -o $build_result -c $config_path
cd $build_result 


# Check for config file before run PSI & dPSI
echo 'START MAJIQ PSI & DPSI'
for index in "${!array[@]}"
do
    treatment=${array[index]}
    majiq_files_treatment=$treatment*.majiq
    psi_path=$sample_path/psi_$treatment

    echo "=================================="
    echo 'PSI for '${treatment} 
    majiq psi -o $psi_path -n $treatment $majiq_files_treatment --output-type all

    if (( "${#array[@]}" < 2 ))
    then
        if [ $index != 0 ]
        then
            control=${array[0]}
            majiq_files_control=${control}*.majiq
            echo "----------------------------------"
            echo 'dPSI for '${control} ${treatment} 
            majiq deltapsi -o /home/dhthutrang/BN/MAJIQ/$group_ids/dPSI_result/${control}_${treatment} -n $control $treatment -grp1 $majiq_files_control -grp2 $majiq_files_treatment --output-type all
        fi
    fi
done
