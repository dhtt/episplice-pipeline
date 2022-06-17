import argparse
import sys
import os
import subprocess
from pathlib import Path

sys.path.append("/home/dhthutrang/Krebs/episplice-pipeline/utils")
import fastq_dir_to_samplesheet #noLint


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p", 
        "--chip_seq_path", 
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
        "-ctl", 
        "--control_id", 
        type=str,
        help="Sample control ID"
    )
    parser.add_argument(
        "-trm", 
        "--treatment_id", 
        type=str,
        help="Sample treatment IDs"
    )
    parser.add_argument(
        "-ab", 
        "--antibody", 
        type=str,
        help="ChIP-seq antibody"
    )
    parser.add_argument(
        "-ri",
        "--replicate_index",
        type=int,
        default=1,
        help="After splitting FastQ file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name.",
    )
    parser.add_argument(
        "-gi",
        "--group_name_index",
        type=str,
        default="0",
        help="After splitting FastQ file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name.",
    )
    parser.add_argument(
        "-nco",
        "--chipseq_options",
        type=str,
        default="0",
        help="After splitting FastQ file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name.",
    )
    return parser.parse_args(args)


def execute_workflow(args=None):
    # os.chdir('/home/dhthutrang/Krebs/episplice-pipeline')
    args = parse_args(args)
    print(args)
    
    chip_seq_path = args.chip_seq_path
    read1_extension = args.read1_extension
    read2_extension = args.read2_extension
    genome = args.genome
    control_id = args.control_id
    treatment_id = args.treatment_id.split(',')
    antibody = str(args.antibody)
    replicate_index = str(args.replicate_index)
    group_name_index = str(args.group_name_index)
    chipseq_options = args.chipseq_options

    # #TEST VALUES
    # chip_seq_path = '/home/dhthutrang/Krebs/Messier/mRNA_seq'
    # read1_extension = '.fastq.gz'
    # read2_extension = '.fastq.gz'
    # genome = 'GRCh38'
    # control_id = 'MCF7_ctrl'
    # treatment_id = ['MCF7_e2_', 'MCF7_gc10_']
    # antibody = 2
    # replicate_index = 3
    # group_name_index = '0,1,2'
    
    fastq_path = '/'.join([chip_seq_path, 'fastq_files']) #Files in fastq_path must follow naming format celltype_condition_repx_read.fastq.gz
    chipseq_output = '/'.join([chip_seq_path, 'chipseq_output'])
    samplesheet_path = '/'.join([fastq_path, 'samplesheet.csv'])
    print ("========== Initialized DHM workflow ==========")
    # Run nf-core/chip_seq
    print("Create output folder")
    Path(chipseq_output).parent.mkdir(parents=True, exist_ok=True)
    
    print("Make sample sheet")
    fastq_dir_to_samplesheet.main(args = [fastq_path, samplesheet_path, 'chipseq',
        '--read1_extension', read1_extension, 
        '--read2_extension', read2_extension,
        '--control_id', control_id,
        '--antibody', antibody,
        '--replicate_index', replicate_index,
        '--group_name_index', group_name_index
        ])
        
    print("Run nf-core/chipseq")
    subprocess.call(
        "nextflow run nf-core/chipseq --input %s --outdir %s --genome %s -profile docker --save_reference true"
        %(samplesheet_path, chipseq_output, genome) +
        chipseq_options, 
        shell=True)
    
    print("========== Finished nf-core/chipseq ==========")


if __name__ == "__main__":
    sys.exit(execute_workflow())
