import os
import sys
import argparse
from modules.utils import get_bams_from_list, get_annotation_resources
from modules.fragmentomics import run_fragmentomic_analysis

main_dir = os.path.dirname(os.path.realpath(__file__))

def get_args():

    parser = argparse.ArgumentParser(prog="NanoFrag", 
        description="Analaysis of DNA fragmentation from liquid biopsy using nanopore")

    parser.add_argument("--tumor_list", dest="tumor_list", type=str, required=True, help="Analaysis of DNA fragmentation from liquid biopsy using nanopore")
    parser.add_argument("--normal_list", dest="normal_list", type=str, required=True)
    parser.add_argument("--reference", dest="reference", type=str, required=True)
    parser.add_argument("--output_dir", dest="output_dir", type=str, required=True)
    parser.add_argument("--threads", dest="threads", type=int, required=True)

    args = parser.parse_args()

    return args



if __name__ == "__main__":
    args = get_args()

    tumor_list = args.tumor_list
    normal_list = args.normal_list
    output_dir = args.output_dir
    genome = args.reference
    threads = args.threads

    tumor_bams = get_bams_from_list(tumor_list)
    normal_bams = get_bams_from_list(normal_list)

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    ann_dict = get_annotation_resources(main_dir)

    run_fragmentomic_analysis(normal_bams, ann_dict, genome, output_dir, threads, window_size=5000000)
    # run_fragmentomic_analysis(tumor_bams, ann_dict, genome, output_dir, threads, window_size=5000000)

