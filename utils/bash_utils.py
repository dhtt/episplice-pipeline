import subprocess
from pathlib import Path
from os import listdir
from os.path import isdir, isfile, join

def mkdir_p(path_to_create):
    Path(path_to_create).mkdir(parents=True, exist_ok=True)


def duplicate_file(current_path, new_path):
    subprocess.call(
        "cp %s %s" %(current_path, new_path),
        shell=True
        )

def sort_bed_file(bed_file):
    subprocess.call(
        "sort -k1,1 -k2,2n %s > %s"
        %(bed_file, bed_file + ".sorted"),
        shell=True
        )
    subprocess.call(
        "mv %s %s"
        %(bed_file + ".sorted", bed_file),
        shell=True
        )

def rm_file(file):
    subprocess.call(
        "rm %s" %(file),
        shell=True
        )

def rename_file(file, new_file):
    subprocess.call(
        "mv %s %s" %(file, new_file),
        shell=True
        )


def merge_bed_file(input_dir: str, output_dir: str, bed_files: list, merged_bed_name: str, 
                merge_options="", awk_options="", remove_formats=""):
    # Define input and outname
    merged_bed_file_name = join(output_dir, merged_bed_name + '.bed')
    bed_files = [join(input_dir, s) for s in bed_files]

    # Join bed files by cat
    subprocess.call(
        "cat %s > %s" %(' '.join(bed_files), merged_bed_file_name),
        shell=True
    )
    
    # Sort bed files by bedtools
    sort_bed_file(merged_bed_file_name)
    subprocess.call(
        "bedtools merge %s -i %s > %s"
        %(merge_options, merged_bed_file_name, merged_bed_file_name + ".merged"),
        shell=True
        )
    subprocess.call(
        "awk -v OFS='\t' '{ print %s}' %s > %s"
        %(awk_options, merged_bed_file_name + ".merged", merged_bed_file_name),
        shell=True
        )

    removed_files = ' '.join([output_dir + "/*."+i for i in remove_formats.split(',')])
    print("Removing file: " + removed_files)
    rm_file(removed_files)

def execute_manorm(histone_type, peak_dir, read_dir, control_id, treatment_id, manorm_output_dir):
    samples = [s for s in listdir(peak_dir) if (histone_type in s)]
    control_sample = [s for s in samples if (control_id in s)][0]
    treatment_sample = [s for s in samples if (treatment_id in s)][0]
    peak_1 = join(peak_dir, control_sample)
    peak_2 = join(peak_dir, treatment_sample)
    read_1 = join(read_dir, control_sample)
    read_2 = join(read_dir, treatment_sample)
    output_folder = "_".join([histone_type, control_id.split('.')[0], treatment_id.split('.')[0]])

    try:
        manorm_script = "manorm --p1 " + peak_1 + " --p2 " + peak_2 + " --r1 " + read_1 + " --r2 " + read_2 + " -o " + join(manorm_output_dir, output_folder)
        subprocess.call(manorm_script, shell=True)
    except:
        print("An exception occured while running MAnorm")
        pass

def intersect_bed_file(refgen_flank_path: str, bed_file: str, intersect_options: str):
    subprocess.call(
        "bedtools intersect -a %s -b %s %s > %s"
        %(refgen_flank_path, bed_file, intersect_options, bed_file + ".annot"),
        shell=True
        )
    rename_file(bed_file + ".annot", bed_file)

def collapse_bed_file(bed_file: str, collapse_options: str):
    subprocess.call(
        "bedtools groupby -i %s %s > %s"
        %(bed_file, collapse_options, bed_file + ".grouped"),
        shell=True
        )
    rename_file(bed_file + ".grouped", bed_file)
