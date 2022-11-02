import argparse
import sys
import os
import subprocess
from pathlib import Path
from utils.bash_utils import *
import pandas as pd

sys.path.append("/home/dhthutrang/Krebs/episplice-pipeline/utils")
import fastq_dir_to_samplesheet #noLint


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p", 
        "--mrna_seq_path", 
        type=str,
        help="Path to mRNA-seq folder."
    )
    parser.add_argument(
        "-r1", 
        "--read1_extension", 
        type=str,
        help="File extension for read 1."
    )
    parser.add_argument(
        "-r2", 
        "--read2_extension", 
        type=str,
        help="File extension for read 2."
    )
    parser.add_argument(
        "-G", 
        "--genome", 
        type=str,
        help="Reference genome version"
    )
    parser.add_argument(
        "-gtf", 
        "--refgen_path", 
        type=str,
        help="Reference genome gtf file"
    )
    parser.add_argument(
        "-ctl", 
        "--control_id", 
        type=str,
        help="Sample control ID"
    )
    parser.add_argument(
        "-trm", 
        "--treatment_ids", 
        type=str,
        help="Sample treatment IDs"
    )
    parser.add_argument(
        "-nco",
        "--mrnaseq_options",
        type=str,
        default="0",
        help="",
    )
    return parser.parse_args(args)


def execute_workflow(args=None):
    # os.chdir('/home/dhthutrang/Krebs/episplice-pipeline')
    args = parse_args(args)
    print(args)
    
    mrna_seq_path = args.mrna_seq_path
    read1_extension = args.read1_extension
    read2_extension = args.read2_extension
    genome = args.genome
    refgen_path = args.refgen_path
    control_id = args.control_id
    treatment_ids = args.treatment_ids.split(',')
    mrnaseq_options = args.mrnaseq_options

    # #TEST VALUES
    # mrna_seq_path = '/home/dhthutrang/Krebs/Messier/mRNA_seq'
    # read1_extension = '.fastq.gz'
    # read2_extension = '.fastq.gz'
    # genome = 'GRCh38'
    # control_id = 'MCF7_ctrl'
    # treatment_ids = ['MCF7_e2_', 'MCF7_gc10_']
    # refgen_path = "$ENCODE_REFGEN/reference_genome.2021_.corrected.gtf"

    fastq_path = '/'.join([mrna_seq_path, 'fastq_files']) #Files in fastq_path must follow naming format celltype_condition_repx_read.fastq.gz
    rnaseq_output = '/'.join([mrna_seq_path, 'rnaseq_output'])
    samplesheet_path = '/'.join([fastq_path, 'samplesheet.csv'])
    sam_path = '/'.join([mrna_seq_path, 'SAM_files']) # File name in count_path must be celltype_condition_repx.txt
    count_path = '/'.join([mrna_seq_path, 'count_files']) # File name in count_path must be celltype_condition_repx.txt
    DEXSeq_output_path = '/'.join([mrna_seq_path, 'DEXSeq_output'])
    dexseq_output_annotated = '/'.join([mrna_seq_path, 'dexseq_output_annotated'])
    
    print ("========== Initialized DEU workflow ==========")
    # # Run nf-core/rna_seq
    # print("Create output folder")
    # Path(rnaseq_output).parent.mkdir(parents=True, exist_ok=True)

    # print("Make sample sheet")
    # fastq_dir_to_samplesheet.main(args = [fastq_path, samplesheet_path, 'rnaseq',
    #     '--read1_extension', read1_extension, 
    #     '--read2_extension', read2_extension])
        
    # print("Run nf-core/rnaseq")
    # subprocess.call(
    #     "nextflow run nf-core/rnaseq --input %s --outdir %s --genome %s -profile docker --save_reference true"
    #     %(samplesheet_path, rnaseq_output, genome) +
    #     " " + mrnaseq_options,
    #     shell=True)
    
    # print("========== Finished nf-core/rnaseq ==========")

    # # Generate exon counts
    # print("Create SAM files folder")
    # Path(sam_path).parent.mkdir(parents=True, exist_ok=True)

    # print("Convert BAM to SAM")
    # subprocess.call(
    #     "utils/bam2sam.sh -i %s/star_salmon -o %s" 
    #     %(rnaseq_output, sam_path),
    #     shell=True)

    print("Generate exon count")
    Path(count_path).mkdir(parents=True, exist_ok=True)
    subprocess.call(
        "DEU_scripts/generate_exon_count.sh -i %s -o %s -g %s" 
        %(sam_path, count_path, refgen_path),
        shell=True)

    print("========== Generated exon counts ==========")

    # DEXSeq run
    print("Start DEXseq analysis")
    for extension in ["", "csv", "html", "r_data", "plot"]:
        Path('/'.join([DEXSeq_output_path, extension])).mkdir(parents=True, exist_ok=True)
    
    if control_id != "NULL":
        for trm in treatment_ids:
            subprocess.call(
                "Rscript DEU_scripts/DEXSeq_analysis.R -i %s -o %s -a %s -b %s -G %s -n 8" 
                %(count_path, DEXSeq_output_path, control_id, trm, refgen_path),
                shell=True)
    else:
        for trm1 in treatment_ids:
            for trm2 in treatment_ids:
                if trm1 != trm2:
                    subprocess.call(
                    "Rscript DEU_scripts/DEXSeq_analysis.R -i %s -o %s -a %s -b %s -G %s -n 8" 
                    %(count_path, DEXSeq_output_path, trm1, trm2, refgen_path),
                    shell=True) #TODO: parallel run
                    
    # mkdir_p(dexseq_output_annotated)
    # for dirpath, dirs, files in os.walk('/'.join([DEXSeq_output_path, 'csv'])):	 
    #     for file in files:
    #         if file.endswith(".csv"):
    #             file_path = dirpath + "/" + file
    #             annotated_file_path = join(dexseq_output_annotated, '_'.join(file.split('_')[:-1]) + ".bed")
    #             df = pd.read_csv(file_path, sep='\t', header=0)
    #             df.iloc[:, lambda df:[0, 1, 2, 4, 6, 8, 9]].to_csv(annotated_file_path, sep='\t', header=False, index=False)

    #             intersect_bed_file(
    #                 refgen_flank_path = refgen_flank_path, 
    #                 bed_file = annotated_file_path, 
    #                 intersect_options = "-wo -loj -bed"
    #                 )
    #             collapse_bed_file(
    #                 bed_file=annotated_file_path,
    #                 collapse_options="-g 1-9 -c 13,14 -o max"
    #                 )
    print("========== Finished DEXseq analysis ==========")


if __name__ == "__main__":
    sys.exit(execute_workflow())
