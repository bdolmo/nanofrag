import os
import sys
import argparse
from modules.utils import get_bams_from_list, set_annotation_resources, set_sample_configuration
from modules.fragmentomics import run_fragmentomic_analysis
from modules.nucleosome import windowed_protection_scores

main_dir = os.path.dirname(os.path.realpath(__file__))

def get_args():

    parser = argparse.ArgumentParser(prog="NanoFrag", 
        description="Analaysis of DNA fragmentation from liquid biopsy using nanopore")

    parser.add_argument("--tumor_list", dest="tumor_list", type=str, required=True, 
        help="Analaysis of DNA fragmentation from liquid biopsy using nanopore")
    parser.add_argument("--normal_list", dest="normal_list", type=str, required=True)
    parser.add_argument("--reference", dest="reference", type=str, required=True)
    parser.add_argument("--output_dir", dest="output_dir", type=str, required=True)
    parser.add_argument("--threads", dest="threads", type=int, required=True)
    parser.add_argument("--skip_fragmentation", dest="skip_fragmentation", 
        help="Skip fragment size histogram calculation", action="store_true")
    parser.add_argument("--skip_wps", dest="skip_wps", 
        help="Skip windowed protection score calculation", action="store_true")

    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = get_args()

    tumor_list = args.tumor_list
    normal_list = args.normal_list
    output_dir = args.output_dir
    genome = args.reference
    threads = args.threads
    skip_fragmentation = args.skip_fragmentation
    skip_wps = args.skip_wps

    tumor_bams = get_bams_from_list(tumor_list)
    normal_bams = get_bams_from_list(normal_list)

    ann_dict = set_annotation_resources(main_dir)
    sample_list =  set_sample_configuration(tumor_bams, normal_bams)

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    if skip_fragmentation == False:
        run_fragmentomic_analysis(sample_list, ann_dict, genome, output_dir, threads, window_size=5000000)
    
    if skip_wps == False:
        windowed_protection_scores(sample_list, ann_dict, output_dir)
    

